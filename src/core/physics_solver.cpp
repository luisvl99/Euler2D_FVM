#include "physics_solver.h"
#include "physics_eulerphysics.h"
#include "mesh_boundarytype.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>

// ===========================================================================
//  Constructor
// ===========================================================================

EulerSolver::EulerSolver(StructuredMesh* mesh, const FlowParameters& params)
    : m_mesh(mesh)
    , m_params(params)
{
    if(!m_mesh)
        throw std::runtime_error("EulerSolver: null mesh pointer.");

    buildBoundaryConditions();
}

// ===========================================================================
//  Configuration setters
// ===========================================================================

void EulerSolver::setFlowParameters(const FlowParameters& params)
{
    m_params = params;
}

void EulerSolver::setTimeScheme(TimeScheme scheme)
{
    m_scheme = scheme;
}

void EulerSolver::setCFL(double cfl)
{
    // This CFL is about twice the classical Courant number, so the
    // linear stability limit is about 2
    if(cfl <= 0.0 || cfl > 2.0)
        std::cerr << "[EulerSolver] WARNING: CFL=" << cfl
                  << " is outside the stable range (0, 2].\n";
    m_CFL = cfl;
}

// ===========================================================================
//  step()
//
//  One complete time step.  Returns the dt used.
//
//  History bookkeeping
//  -------------------
//  After computing dt (step 2), we:
//    - push dt to m_dtHistory   → parallel with m_residualHistory
//    - add dt to m_simTime      → accumulated physical time
//  Both are cleared by clearHistory().
// ===========================================================================

double EulerSolver::step()
{
    // 1. Residuals
    computeResiduals();

    // 2. CFL time step — clamped so we land exactly on m_tEnd
    double dt = computeTimeStep();
    dt        = std::min(dt, m_tEnd - m_simTime);
    m_lastDt  = dt;

    // 3. Record history (residual norm + dt + accumulated time)
    double resNorm = 0.0;
    {
        double sum = 0.0;
        for(const Cell& c : m_mesh->cells)
        {
            const EulerState& R = c.R;
            sum += R.rho*R.rho + R.rho_u*R.rho_u
                   + R.rho_v*R.rho_v + R.rho_E*R.rho_E;
        }
        resNorm = std::sqrt(sum / static_cast<double>(m_mesh->cells.size()));
    }
    m_residualHistory.push_back(resNorm);
    m_dtHistory.push_back(dt);    // ← parallel array: same index as residual
    m_simTime += dt;              // ← accumulated physical time
    ++m_totalIter;

    // 4. Advance, then stop if any cell became non-physical
    if(m_scheme == TimeScheme::RK2)
        applyRK2(dt);
    else
        applyForwardEuler(dt);

    checkStates();

    return dt;
}

// ===========================================================================
//  Diagnostics
// ===========================================================================

void EulerSolver::clearHistory()
{
    m_residualHistory.clear();
    m_dtHistory.clear();
    m_simTime   = 0.0;
    m_totalIter = 0;
    m_lastDt    = 0.0;
}

void EulerSolver::setReconstructionScheme(ReconstructionScheme r)
{
    m_reconstruction = r;
}
void EulerSolver::setFluxScheme(FluxScheme f)
{
    m_fluxScheme = f;
}

// ---------------------------------------------------------------------------
//  minmod limiter applied component-wise to an EulerState
//
//  minmod(a, b):
//    - if a and b have opposite signs → 0  (be conservative near extrema)
//    - otherwise                      → the one with the smaller magnitude
//
//  Applied to each conserved variable independently.
// ---------------------------------------------------------------------------
static EulerState minmod(const EulerState& a, const EulerState& b)
{
    auto mm = [](double x, double y) -> double {
        if(x * y <= 0.0) return 0.0;
        return (std::abs(x) < std::abs(y)) ? x : y;
    };
    return { mm(a.rho,   b.rho),
            mm(a.rho_u, b.rho_u),
            mm(a.rho_v, b.rho_v),
            mm(a.rho_E, b.rho_E) };
}

// ---------------------------------------------------------------------------
//  MUSCL reconstruction at one face
//
//  Given the four-cell stencil  [LL] | [L] | [R] | [RR]
//  and flags indicating which extended neighbours exist (boundary cells
//  have no LL or RR), returns the reconstructed left/right face states.
//
//  On a uniform mesh the cell spacing cancels, so slopes are in units of
//  state difference (not state difference / dx).  The factor 0.5 comes
//  from the half-cell extrapolation to the face.
// ---------------------------------------------------------------------------
static void musclReconstruct(const EulerState& ULL, const EulerState& UL,
                             const EulerState& UR,  const EulerState& URR,
                             bool hasLL, bool hasRR,
                             EulerState& UL_face, EulerState& UR_face)
{
    // Left cell slope:  difference across [LL→L] and [L→R]
    EulerState slopeL{};
    if(hasLL)
        slopeL = minmod(UL - ULL, UR - UL);
    // else: boundary → zero slope (piecewise constant fallback)

    // Right cell slope:  difference across [L→R] and [R→RR]
    EulerState slopeR{};
    if(hasRR)
        slopeR = minmod(UR - UL, URR - UR);

    UL_face = UL + slopeL * 0.5;   // extrapolate LEFT  cell to face
    UR_face = UR - slopeR * 0.5;   // extrapolate RIGHT cell to face
}



// ===========================================================================
//  computeResiduals
// ===========================================================================

void EulerSolver::computeResiduals()
{
    for(Cell& c : m_mesh->cells)
        c.R = EulerState{};

    for(const Face& face : m_mesh->faces)
    {
        const int L = face.leftCell;
        const int R = face.rightCell;

        if(face.boundaryType == BoundaryType::Internal)
        {
            EulerState UL_face, UR_face;

            if(m_reconstruction == ReconstructionScheme::MUSCL)
            {
                const Cell& cellL = m_mesh->cells[L];
                const Cell& cellR = m_mesh->cells[R];

                // Determine face orientation and look up the extended stencil.
                // The structured mesh stores .i and .j in every Cell, so we can
                // directly address neighbours without a separate adjacency list.
                EulerState ULL{}, URR{};
                bool hasLL = false, hasRR = false;

                if(cellR.i == cellL.i + 1)          // vertical face  (x-sweep)
                {
                    if(cellL.i - 1 >= 0) {
                        ULL   = m_mesh->cells[m_mesh->cellIndex(cellL.i - 1, cellL.j)].U;
                        hasLL = true;
                    }
                    if(cellR.i + 1 < m_mesh->Nx) {
                        URR   = m_mesh->cells[m_mesh->cellIndex(cellR.i + 1, cellR.j)].U;
                        hasRR = true;
                    }
                }
                else                                 // horizontal face (y-sweep)
                {
                    if(cellL.j - 1 >= 0) {
                        ULL   = m_mesh->cells[m_mesh->cellIndex(cellL.i, cellL.j - 1)].U;
                        hasLL = true;
                    }
                    if(cellR.j + 1 < m_mesh->Ny) {
                        URR   = m_mesh->cells[m_mesh->cellIndex(cellR.i, cellR.j + 1)].U;
                        hasRR = true;
                    }
                }

                musclReconstruct(ULL, cellL.U, cellR.U, URR,
                                 hasLL, hasRR,
                                 UL_face, UR_face);
            }
            else
            {
                UL_face = m_mesh->cells[L].U;
                UR_face = m_mesh->cells[R].U;
            }

            EulerState flux = (m_fluxScheme == FluxScheme::HLLC
                                   ? EulerPhysics::hllcFlux(UL_face, UR_face, face.nx, face.ny)
                                   : EulerPhysics::rusanovFlux(UL_face, UR_face, face.nx, face.ny))
                              * face.length;
            m_mesh->cells[L].R -= flux;
            m_mesh->cells[R].R += flux;
        }
        else
        {
            const int interiorIdx = (L >= 0) ? L : R;
            const EulerState& U_real = m_mesh->cells[interiorIdx].U;

            auto it = m_bc.find(static_cast<int>(face.boundaryType));
            if(it == m_bc.end())
            {
                std::cerr << "[EulerSolver] WARNING: no BC handler for type "
                          << static_cast<int>(face.boundaryType) << "\n";
                continue;
            }

            EulerState U_ghost = it->second->ghostState(U_real, face, m_params);
            EulerState flux = (m_fluxScheme == FluxScheme::HLLC
                                   ? EulerPhysics::hllcFlux(U_real, U_ghost, face.nx, face.ny)
                                   : EulerPhysics::rusanovFlux(U_real, U_ghost, face.nx, face.ny))
                              * face.length;
            m_mesh->cells[interiorIdx].R -= flux;
        }
    }
}

// ===========================================================================
//  computeTimeStep
// ===========================================================================

double EulerSolver::computeTimeStep() const
{
    double dt_min = std::numeric_limits<double>::max();

    for(const Cell& cell : m_mesh->cells)
    {
        double sum_lambda = 0.0;

        for(int fIdx : cell.faces)
        {
            const Face& face = m_mesh->faces[fIdx];
            double lam = EulerPhysics::maxWaveSpeed(cell.U, face.nx, face.ny);
            sum_lambda += lam * face.length;
        }

        if(sum_lambda > 1e-30)
        {
            double dt_cell = m_CFL * cell.area / sum_lambda;
            dt_min = std::min(dt_min, dt_cell);
        }
    }

    if(dt_min == std::numeric_limits<double>::max())
    {
        std::cerr << "[EulerSolver] WARNING: could not compute finite dt.\n";
        dt_min = 1e-6;
    }

    return dt_min;
}

// ===========================================================================
//  applyForwardEuler
// ===========================================================================

void EulerSolver::applyForwardEuler(double dt)
{
    for(Cell& cell : m_mesh->cells)
        cell.U += cell.R * (dt / cell.area);
}

// ===========================================================================
//  applyRK2  (SSP-RK2 / Heun)
// ===========================================================================

void EulerSolver::applyRK2(double dt)
{
    // SSP-RK2 (Heun / Shu-Osher):
    //
    //   Stage 1:  U*    = U^n  +  dt * R^n / Ω
    //   Stage 2:  U^n+1 = 0.5 * (U^n  +  U*  +  dt * R* / Ω)
    //
    // cell.R holds R^n on entry (written by the computeResiduals() call in
    // step(), which also recorded the residual-history entry).
    // We must NOT let the stage-2 computeResiduals() permanently overwrite
    // cell.R, because the viewer reads cell.R to display the "Residual |R|"
    // field — it must stay at R^n so it is consistent with the history plot.

    const std::size_t N = m_mesh->cells.size();
    m_U0.resize(N);
    m_R0.resize(N);

    // --- save U^n and R^n; advance to U* ---
    for(std::size_t i = 0; i < N; ++i)
    {
        Cell& cell = m_mesh->cells[i];
        m_U0[i] = cell.U;
        m_R0[i] = cell.R;                             // save R^n
        cell.U  = m_U0[i] + cell.R * (dt / cell.area);
    }

    checkStates();        // U* must be physical before R* is evaluated

    // --- R* = residual at U* ---
    computeResiduals();   // writes cell.R = R*

    // --- complete the SSP average; restore R^n into cell.R ---
    for(std::size_t i = 0; i < N; ++i)
    {
        Cell& cell  = m_mesh->cells[i];
        EulerState U_ss = cell.U + cell.R * (dt / cell.area);
        cell.U = (m_U0[i] + U_ss) * 0.5;
        cell.R = m_R0[i];                             // restore R^n
    }
}

// ===========================================================================
//  buildBoundaryConditions
// ===========================================================================

void EulerSolver::buildBoundaryConditions()
{
    for(const Face& face : m_mesh->faces)
    {
        const int key = static_cast<int>(face.boundaryType);

        if(face.boundaryType == BoundaryType::Internal) continue;

        if(m_bc.count(key) == 0)
            m_bc[key] = BoundaryConditionFactory::create(face.boundaryType);
    }
}

// ===========================================================================
//  checkStates
//
//  Throws std::runtime_error at the first cell whose state is not
//  physically admissible (rho > 0 and p > 0, both finite).  Called after
//  every update, so a blow-up stops the run at the step where it starts
//  instead of spreading NaN through the field.
// ===========================================================================

void EulerSolver::checkStates() const
{
    for(const Cell& c : m_mesh->cells)
    {
        const double p = EulerPhysics::pressure(c.U);
        if(c.U.rho > 0.0 && p > 0.0 && std::isfinite(c.U.rho) && std::isfinite(p))
            continue;

        std::ostringstream msg;
        msg << "Non-physical state in cell (i=" << c.i << ", j=" << c.j
            << ") at step " << m_totalIter << ", t = " << m_simTime
            << ": rho = " << c.U.rho << ", p = " << p
            << ". Try a smaller CFL number.";
        throw std::runtime_error(msg.str());
    }
}
