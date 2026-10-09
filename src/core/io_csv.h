#pragma once

#include <iosfwd>
#include <string>
#include <vector>

#include "mesh_structured.h"
#include "physics_eulerstate.h"
#include "physics_flowparameters.h"
#include "physics_solver.h"
#include "physics_testcases.h"

// ---------------------------------------------------------------------------
// CSV output shared by the GUI and the command-line runner
//
//  Every file starts with '#' metadata lines describing the run, so pandas
//  reads it directly with  pd.read_csv(path, comment="#").
// ---------------------------------------------------------------------------

// Everything the metadata header records about a run
struct RunInfo
{
    std::string caseName;                  // testCaseId(), e.g. "sod"
    bool        subsonic = false;
    double      tEnd     = 0.0;

    int    Nx = 0, Ny = 0;
    double Lx = 0.0, Ly = 0.0;

    double rho_inf = 0.0, u_inf = 0.0, v_inf = 0.0, p_inf = 0.0, p_back = 0.0;

    double      CFL = 0.0;
    std::string timeScheme;                // "ForwardEuler" or "RK2"
    std::string flux;                      // "Rusanov" or "HLLC"
    std::string reconstruction;            // "PC" or "MUSCL"

    int maxIter       = 0;
    int snapshotEvery = 0;                 // 0 = no periodic snapshots
};

// Fills RunInfo from the objects actually used by the run (maxIter and
// snapshotEvery are left for the caller)
RunInfo makeRunInfo(TestCase tc, bool subsonic, const StructuredMesh& mesh,
                    const FlowParameters& fp, const EulerSolver& solver);

// Columns: iter, dt, sim_time, residual_L2.  Row k is the residual of the
// state after k steps (at time sim_time) and the dt of the step taken from it.
void writeResidualsCsv(std::ostream& out, const RunInfo& info,
                       const std::vector<double>& residuals,
                       const std::vector<double>& dts);

// Columns: i, j, cx, cy, rho, u, v, p, mach, rho_E.  U[k] belongs to mesh.cells[k].
void writeSnapshotCsv(std::ostream& out, const RunInfo& info,
                      const StructuredMesh& mesh, const std::vector<EulerState>& U,
                      int iter, double simTime);

// Columns: i, x, rho, u, v, p, mach for the cells of row j (0 <= j < Ny)
void writeLineProbeCsv(std::ostream& out, const RunInfo& info,
                       const StructuredMesh& mesh, const std::vector<EulerState>& U,
                       int j, int iter, double simTime);
