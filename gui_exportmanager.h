#pragma once

#include <QString>
#include <vector>

#include "physics_eulerstate.h"
#include "mesh_structured.h"

// ---------------------------------------------------------------------------
// SimParameters
//
//  Single record of every user-controlled parameter active when a run was
//  started.  Captured once in MainWindow::on_btnRunSolver_clicked() and
//  attached to every exported file as metadata comment lines.
//
//  Python can read these with:
//    pd.read_csv("file.csv", comment='#')
// ---------------------------------------------------------------------------
struct SimParameters
{
    // Mesh
    int    Nx = 0,  Ny = 0;
    double Lx = 0., Ly = 0.;

    // Free-stream / inlet
    double rho_inf = 0., u_inf = 0., v_inf = 0., p_inf = 0.;

    // Solver
    double  CFL          = 0.5;
    QString scheme;           // "RK2" or "ForwardEuler"
    int     maxIter      = 0;
    int     snapshotEvery = 0;

    QString fluxScheme;       // "Rusanov" or "HLLC"
    QString reconstruction;   // "PC" or "MUSCL"
};

// ---------------------------------------------------------------------------
// CellSnapshot
//
//  A full copy of every cell's conserved state at one point in time.
//  Stored by MainWindow every snapshotEvery iterations.
//  Passed to ExportManager without touching the live mesh.
//
//  cells[k]  corresponds to  mesh->cells[k]  (same flat ordering).
// ---------------------------------------------------------------------------
struct CellSnapshot
{
    int    iter    = 0;
    double simTime = 0.0;
    std::vector<EulerState> cells;   // one EulerState per mesh cell
};

// ---------------------------------------------------------------------------
// ExportManager
//
//  Pure-static utility class.  All write functions return true on success.
//  On failure they return false and, if errorOut != nullptr, set it to a
//  human-readable message suitable for a QMessageBox.
//
//  All files share the same metadata header format:
//
//    # ==================================================
//    # Euler2D_FVM — <title>
//    # ==================================================
//    # Exported    : <ISO datetime>
//    # Mesh        : Nx=…  Ny=…  Lx=…  Ly=…
//    # Flow (inf)  : rho=…  u=…  v=…  p=…
//    # Solver      : CFL=…  Scheme=…  MaxIter=…
//    # Snapshot    : iter=…  t=…  [snapshotEvery=…]
//
//  Python usage:
//    import pandas as pd
//    df = pd.read_csv("file.csv", comment='#')
// ---------------------------------------------------------------------------
class ExportManager
{
public:

    // -----------------------------------------------------------------------
    // writeResiduals
    //
    //  Columns: iter, dt, sim_time, residual_L2
    //  One row per solver step (full history).
    // -----------------------------------------------------------------------
    static bool writeResiduals(const QString&             path,
                               const std::vector<double>& residuals,
                               const std::vector<double>& dtHistory,
                               const SimParameters&       params,
                               QString*                   errorOut = nullptr);

    // -----------------------------------------------------------------------
    // writeSnapshot
    //
    //  Columns: i, j, cx, cy, rho, u, v, p, mach, rho_E
    //  One row per cell, for the conserved state stored in 'snap'.
    // -----------------------------------------------------------------------
    static bool writeSnapshot(const QString&        path,
                              const StructuredMesh* mesh,
                              const CellSnapshot&   snap,
                              const SimParameters&  params,
                              QString*              errorOut = nullptr);

    // -----------------------------------------------------------------------
    // writeLineProbe
    //
    //  1-D slice at fixed j row.
    //  Columns: i, x, rho, u, v, p, mach
    //  One row per i index (all cells at fixed j, ordered i=0…Nx-1).
    // -----------------------------------------------------------------------
    static bool writeLineProbe(const QString&        path,
                               const StructuredMesh* mesh,
                               const CellSnapshot&   snap,
                               int                   jRow,
                               const SimParameters&  params,
                               QString*              errorOut = nullptr);

private:

    // Build the standard metadata comment block (shared by all write functions)
    static QString makeHeader(const SimParameters& params,
                              const QString&       title,
                              int                  iter,
                              double               simTime);
};
