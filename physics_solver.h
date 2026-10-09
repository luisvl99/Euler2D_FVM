#pragma once

#include <memory>
#include <unordered_map>
#include <vector>

#include "mesh_structured.h"
#include "physics_boundarycondition.h"
#include "physics_eulerstate.h"
#include "physics_flowparameters.h"

// ---------------------------------------------------------------------------
// EulerSolver
//
//  Explicit finite-volume solver for the 2-D Euler equations on a
//  structured quadrilateral mesh.
//
//  Numerical method
//  ----------------
//  * Numerical flux         : Rusanov (Local Lax-Friedrichs)  or  HLLC
//  * Reconstruction         : piecewise constant  or  MUSCL + MinMod
//  * Time integration       : Forward Euler  or  SSP-RK2 (Heun)
//  * Boundary conditions    : ghost-cell approach via BoundaryCondition hierarchy
//  * Time step              : global CFL condition  dt = CFL * min_i(Ω_i / Σλ_f)
//                             (sum over all faces, so this CFL is about twice
//                             the classical Courant number)
//
//  Ownership
//  ---------
//  The solver holds a NON-OWNING pointer to the mesh.  The mesh must outlive
//  the solver.  FlowParameters is copied by value so the solver is
//  self-contained once constructed.
// ---------------------------------------------------------------------------

class EulerSolver
{
public:

    // -----------------------------------------------------------------------
    // Time integration schemes
    // -----------------------------------------------------------------------
    enum class TimeScheme
    {
        ForwardEuler,
        RK2
    };

    enum class ReconstructionScheme
    {
        PiecewiseConstant,   // cell average up to the face — first order
        MUSCL                // piecewise linear + MinMod limiter — second order
                             // in smooth regions (with SSP-RK2)
    };

    enum class FluxScheme
    {
        Rusanov,
        HLLC
    };

    // -----------------------------------------------------------------------
    // Constructor
    // -----------------------------------------------------------------------
    EulerSolver(StructuredMesh* mesh, const FlowParameters& params);

    // -----------------------------------------------------------------------
    // Configuration
    // -----------------------------------------------------------------------
    void setFlowParameters(const FlowParameters& params);
    void setTimeScheme(TimeScheme scheme);
    void setCFL(double cfl);
    void setReconstructionScheme(ReconstructionScheme r);
    void setFluxScheme(FluxScheme f);

    // Stop time: step() clamps dt so simTime never exceeds this value.
    // Default: infinity (run forever).
    void   setTEnd(double t) { m_tEnd = t; }
    double tEnd()      const { return m_tEnd; }

    double     CFL()    const { return m_CFL;    }
    TimeScheme scheme() const { return m_scheme; }

    // -----------------------------------------------------------------------
    // Time advancement
    // -----------------------------------------------------------------------
    // One step; returns the dt used.  Throws std::runtime_error if a cell
    // becomes non-physical (rho or p not positive / not finite).
    double step();

    // -----------------------------------------------------------------------
    // Diagnostics
    // -----------------------------------------------------------------------
    void   clearHistory();       // clears residuals, dtHistory, resets simTime

    // Per-step L2 residual norm (one entry per step())
    const std::vector<double>& residualHistory() const { return m_residualHistory; }

    // Per-step dt (same length as residualHistory — parallel arrays)
    const std::vector<double>& dtHistory() const { return m_dtHistory; }

    int    totalIterations() const { return m_totalIter; }
    double lastDt()          const { return m_lastDt;   }

    // Accumulated physical time:  simTime = Σ dt  over all steps
    double simTime() const { return m_simTime; }

private:

    // -----------------------------------------------------------------------
    // Numerical kernels
    // -----------------------------------------------------------------------
    void   computeResiduals();
    double computeTimeStep() const;
    void   applyForwardEuler(double dt);
    void   applyRK2(double dt);

    // -----------------------------------------------------------------------
    // Internal helpers
    // -----------------------------------------------------------------------
    void buildBoundaryConditions();
    void checkStates() const;   // throws at the first non-physical cell

    // -----------------------------------------------------------------------
    // Data
    // -----------------------------------------------------------------------
    StructuredMesh* m_mesh   = nullptr;
    FlowParameters  m_params;
    TimeScheme      m_scheme = TimeScheme::RK2;
    double          m_CFL    = 0.5;
    double          m_tEnd   = 1e30;  // physical stop time (clamped in step())

    ReconstructionScheme m_reconstruction = ReconstructionScheme::PiecewiseConstant;
    FluxScheme m_fluxScheme = FluxScheme::Rusanov;

    std::unordered_map<int, std::shared_ptr<BoundaryCondition>> m_bc;
    std::vector<EulerState> m_U0;   // RK2 stage-1 conserved state  (U^n)
    std::vector<EulerState> m_R0;   // RK2 stage-1 residual         (R^n)

    // History (parallel arrays — same index = same step)
    std::vector<double> m_residualHistory;
    std::vector<double> m_dtHistory;

    int    m_totalIter = 0;
    double m_lastDt    = 0.0;
    double m_simTime   = 0.0;   // Σ dt — accumulated physical time
};
