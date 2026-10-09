#include "io_csv.h"
#include "physics_eulerphysics.h"

#include <cmath>
#include <ctime>
#include <iomanip>
#include <ostream>
#include <sstream>
#include <stdexcept>

namespace {

const char* toString(EulerSolver::TimeScheme s)
{
    return s == EulerSolver::TimeScheme::RK2 ? "RK2" : "ForwardEuler";
}

const char* toString(EulerSolver::FluxScheme f)
{
    return f == EulerSolver::FluxScheme::HLLC ? "HLLC" : "Rusanov";
}

const char* toString(EulerSolver::ReconstructionScheme r)
{
    return r == EulerSolver::ReconstructionScheme::MUSCL ? "MUSCL" : "PC";
}

// '#' metadata block that starts every file
void writeHeader(std::ostream& out, const RunInfo& r, const std::string& title,
                 int iter, double simTime)
{
    const char* sep = "# ============================================================\n";
    const std::time_t now = std::time(nullptr);

    out << std::setprecision(10)
        << sep
        << "# Euler2D_FVM - " << title << "\n"
        << sep
        << "# Exported    : " << std::put_time(std::localtime(&now), "%Y-%m-%d %H:%M:%S") << "\n"
        << "#\n"
        << "# --- Case -------------------------------------------------------\n"
        << "# case=" << r.caseName << "  t_end=" << r.tEnd
        << "  subsonic_bcs=" << (r.subsonic ? 1 : 0) << "\n"
        << "#\n"
        << "# --- Mesh -------------------------------------------------------\n"
        << "# Nx=" << r.Nx << "  Ny=" << r.Ny << "  Lx=" << r.Lx << "  Ly=" << r.Ly << "\n"
        << "#\n"
        << "# --- Flow (free-stream / inlet) ---------------------------------\n"
        << "# rho=" << r.rho_inf << "  u=" << r.u_inf << "  v=" << r.v_inf
        << "  p=" << r.p_inf << "  p_back=" << r.p_back
        << "  gamma=" << EulerPhysics::gamma << "\n"
        << "#\n"
        << "# --- Solver -----------------------------------------------------\n"
        << "# CFL=" << r.CFL << "  Scheme=" << r.timeScheme << "  Flux=" << r.flux
        << "  Reconstruction=" << r.reconstruction << "  MaxIter=" << r.maxIter << "\n"
        << "#\n"
        << "# --- Snapshot ---------------------------------------------------\n"
        << "# iter=" << iter << "  sim_time=" << simTime;
    if(r.snapshotEvery > 0)
        out << "  snapshot_every=" << r.snapshotEvery;
    out << "\n" << sep;
}

struct Primitives { double rho, u, v, p, mach; };

Primitives primitives(const EulerState& U)
{
    const double rho = U.rho;
    const double u   = (rho > 1e-30) ? U.rho_u / rho : 0.0;
    const double v   = (rho > 1e-30) ? U.rho_v / rho : 0.0;
    const double a   = EulerPhysics::soundSpeed(U);
    const double M   = (a > 1e-30) ? std::sqrt(u*u + v*v) / a : 0.0;
    return {rho, u, v, EulerPhysics::pressure(U), M};
}

void requireOneStatePerCell(const StructuredMesh& mesh, const std::vector<EulerState>& U)
{
    if(U.size() != mesh.cells.size())
        throw std::invalid_argument("CSV export: state count does not match the mesh.");
}

} // anonymous namespace

RunInfo makeRunInfo(TestCase tc, bool subsonic, const StructuredMesh& mesh,
                    const FlowParameters& fp, const EulerSolver& solver)
{
    RunInfo r;
    r.caseName = testCaseId(tc);
    r.subsonic = subsonic;
    r.tEnd     = solver.tEnd();

    r.Nx = mesh.Nx;  r.Ny = mesh.Ny;
    r.Lx = mesh.Lx;  r.Ly = mesh.Ly;

    r.rho_inf = fp.rho_inf;  r.u_inf = fp.u_inf;  r.v_inf = fp.v_inf;
    r.p_inf   = fp.p_inf;    r.p_back = fp.p_back;

    r.CFL            = solver.CFL();
    r.timeScheme     = toString(solver.scheme());
    r.flux           = toString(solver.fluxScheme());
    r.reconstruction = toString(solver.reconstruction());
    return r;
}

void writeResidualsCsv(std::ostream& out, const RunInfo& info,
                       const std::vector<double>& residuals,
                       const std::vector<double>& dts)
{
    double tFinal = 0.0;
    for(double dt : dts)
        tFinal += dt;

    writeHeader(out, info, "Residual History", static_cast<int>(residuals.size()), tFinal);
    out << "iter,dt,sim_time,residual_L2\n";

    double t = 0.0;
    for(std::size_t k = 0; k < residuals.size(); ++k)
    {
        const double dt = (k < dts.size()) ? dts[k] : 0.0;
        out << k << ',' << dt << ',' << t << ',' << residuals[k] << '\n';
        t += dt;
    }
}

void writeSnapshotCsv(std::ostream& out, const RunInfo& info,
                      const StructuredMesh& mesh, const std::vector<EulerState>& U,
                      int iter, double simTime)
{
    requireOneStatePerCell(mesh, U);

    writeHeader(out, info, "Cell Snapshot", iter, simTime);
    out << "i,j,cx,cy,rho,u,v,p,mach,rho_E\n";

    for(std::size_t k = 0; k < mesh.cells.size(); ++k)
    {
        const Cell&      c = mesh.cells[k];
        const Primitives w = primitives(U[k]);
        out << c.i  << ',' << c.j  << ',' << c.cx << ',' << c.cy << ','
            << w.rho << ',' << w.u << ',' << w.v  << ',' << w.p  << ','
            << w.mach << ',' << U[k].rho_E << '\n';
    }
}

void writeLineProbeCsv(std::ostream& out, const RunInfo& info,
                       const StructuredMesh& mesh, const std::vector<EulerState>& U,
                       int j, int iter, double simTime)
{
    requireOneStatePerCell(mesh, U);
    if(j < 0 || j >= mesh.Ny)
        throw std::invalid_argument("CSV export: line probe row out of range.");

    std::ostringstream title;
    title << "Line Probe  j=" << j << " (y=" << mesh.cells[mesh.cellIndex(0, j)].cy << ")";

    writeHeader(out, info, title.str(), iter, simTime);
    out << "i,x,rho,u,v,p,mach\n";

    for(int i = 0; i < mesh.Nx; ++i)
    {
        const int        k = mesh.cellIndex(i, j);
        const Primitives w = primitives(U[k]);
        out << i << ',' << mesh.cells[k].cx << ','
            << w.rho << ',' << w.u << ',' << w.v << ',' << w.p << ',' << w.mach << '\n';
    }
}
