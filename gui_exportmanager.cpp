#include "gui_exportmanager.h"
#include "physics_eulerphysics.h"

#include <QFile>
#include <QTextStream>
#include <QDateTime>

#include <cmath>

// ===========================================================================
//  Internal helper: open a file and return a ready QTextStream
//  Returns false (sets errorOut) if the file cannot be opened.
// ===========================================================================
static bool openFile(const QString& path,
                     QFile& file,
                     QTextStream& ts,
                     QString* errorOut)
{
    file.setFileName(path);
    if(!file.open(QIODevice::WriteOnly | QIODevice::Text))
    {
        if(errorOut)
            *errorOut = QString("Cannot open file for writing:\n%1\n\n%2")
                            .arg(path, file.errorString());
        return false;
    }
    ts.setDevice(&file);
    ts.setRealNumberPrecision(10);
    return true;
}

// ===========================================================================
//  makeHeader
//
//  Produces the standard metadata block that starts every exported file.
//  All lines begin with '#' so Python's pd.read_csv(..., comment='#')
//  skips them automatically.
// ===========================================================================
QString ExportManager::makeHeader(const SimParameters& p,
                                  const QString&       title,
                                  int                  iter,
                                  double               simTime)
{
    QString h;
    QTextStream s(&h);
    QString sep = "# ============================================================\n";

    s << sep;
    s << "# Euler2D_FVM — " << title << "\n";
    s << sep;
    s << "# Exported    : "
      << QDateTime::currentDateTime().toString("yyyy-MM-dd hh:mm:ss") << "\n";
    s << "#\n";
    s << "# --- Mesh -------------------------------------------------------\n";
    s << "# Nx=" << p.Nx << "  Ny=" << p.Ny
      << "  Lx=" << p.Lx << "  Ly=" << p.Ly << "\n";
    s << "#\n";
    s << "# --- Flow (free-stream / inlet) ---------------------------------\n";
    s << "# rho=" << p.rho_inf
      << "  u=" << p.u_inf
      << "  v=" << p.v_inf
      << "  p=" << p.p_inf << "\n";
    s << "#\n";
    s << "# --- Solver -----------------------------------------------------\n";
    s << "# CFL=" << p.CFL
      << "  Scheme=" << p.scheme
      << "  Flux=" << p.fluxScheme
      << "  Reconstruction=" << p.reconstruction
      << "  MaxIter=" << p.maxIter << "\n";
    s << "#\n";
    s << "# --- Snapshot ---------------------------------------------------\n";
    s << "# iter="    << iter
      << "  sim_time=" << simTime;
    if(p.snapshotEvery > 0)
        s << "  snapshot_every=" << p.snapshotEvery;
    s << "\n";
    s << sep;

    return h;
}

// ===========================================================================
//  writeResiduals
//
//  Columns: iter, dt, sim_time, residual_L2
//
//  Both residuals and dtHistory must be the same length (they are parallel
//  arrays filled together in EulerSolver::step()).  sim_time is the prefix
//  sum of dt, computed here so it does not need to be stored separately.
// ===========================================================================
bool ExportManager::writeResiduals(const QString&             path,
                                   const std::vector<double>& residuals,
                                   const std::vector<double>& dtHistory,
                                   const SimParameters&       params,
                                   QString*                   errorOut)
{
    if(residuals.empty())
    {
        if(errorOut) *errorOut = "No residual data to export (run the solver first).";
        return false;
    }

    QFile file;
    QTextStream ts;
    if(!openFile(path, file, ts, errorOut)) return false;

    // Metadata header
    int    lastIter = static_cast<int>(residuals.size()) - 1;
    double totalT   = 0.0;
    for(double d : dtHistory) totalT += d;

    ts << makeHeader(params, "Residual History", lastIter, totalT);
    ts << "iter,dt,sim_time,residual_L2\n";

    const int    N  = static_cast<int>(residuals.size());
    const int    Nd = static_cast<int>(dtHistory.size());
    double       t  = 0.0;

    for(int i = 0; i < N; ++i)
    {
        double dt = (i < Nd) ? dtHistory[i] : 0.0;
        t += dt;
        ts << i << ","
           << dt << ","
           << t  << ","
           << residuals[i] << "\n";
    }

    return true;
}

// ===========================================================================
//  writeSnapshot
//
//  Columns: i, j, cx, cy, rho, u, v, p, mach, rho_E
//
//  All primitive variables are derived from the stored conserved state using
//  EulerPhysics, so the file is self-consistent even if the reader does not
//  know the EOS.
// ===========================================================================
bool ExportManager::writeSnapshot(const QString&        path,
                                  const StructuredMesh* mesh,
                                  const CellSnapshot&   snap,
                                  const SimParameters&  params,
                                  QString*              errorOut)
{
    if(!mesh || snap.cells.empty())
    {
        if(errorOut) *errorOut = "No snapshot data available.";
        return false;
    }

    QFile file;
    QTextStream ts;
    if(!openFile(path, file, ts, errorOut)) return false;

    ts << makeHeader(params, "Cell Snapshot", snap.iter, snap.simTime);
    ts << "i,j,cx,cy,rho,u,v,p,mach,rho_E\n";

    const int Nx = mesh->Nx;
    const int Ny = mesh->Ny;

    for(int j = 0; j < Ny; ++j)
    {
        for(int i = 0; i < Nx; ++i)
        {
            int idx = j * Nx + i;   // matches mesh->cellIndex(i,j)

            const Cell&       cell = mesh->cells[idx];
            const EulerState& U    = snap.cells[idx];

            // Primitive variables
            double rho = U.rho;
            double u   = (rho > 1e-30) ? U.rho_u / rho : 0.0;
            double v   = (rho > 1e-30) ? U.rho_v / rho : 0.0;
            double p   = EulerPhysics::pressure(U);
            double a   = EulerPhysics::soundSpeed(U);
            double spd = std::sqrt(u*u + v*v);
            double M   = (a > 1e-30) ? spd / a : 0.0;

            ts << i         << ","
               << j         << ","
               << cell.cx   << ","
               << cell.cy   << ","
               << rho        << ","
               << u          << ","
               << v          << ","
               << p          << ","
               << M          << ","
               << U.rho_E    << "\n";
        }
    }

    return true;
}

// ===========================================================================
//  writeLineProbe
//
//  Columns: i, x, rho, u, v, p, mach
//
//  Extracts a 1-D slice at fixed j=jRow, sweeping i=0…Nx-1.
//  'x' is the cell-centre x-coordinate.  This directly produces a file
//  that can be overlaid with the exact Riemann / Sod solution in Python.
// ===========================================================================
bool ExportManager::writeLineProbe(const QString&        path,
                                   const StructuredMesh* mesh,
                                   const CellSnapshot&   snap,
                                   int                   jRow,
                                   const SimParameters&  params,
                                   QString*              errorOut)
{
    if(!mesh || snap.cells.empty())
    {
        if(errorOut) *errorOut = "No snapshot data available.";
        return false;
    }
    if(jRow < 0 || jRow >= mesh->Ny)
    {
        if(errorOut) *errorOut = QString("j=%1 is out of range [0, %2].")
                            .arg(jRow).arg(mesh->Ny - 1);
        return false;
    }

    QFile file;
    QTextStream ts;
    if(!openFile(path, file, ts, errorOut)) return false;

    // Title encodes the j row so files are self-documenting
    QString title = QString("Line Probe  j=%1 (y=%2)")
                        .arg(jRow)
                        .arg(mesh->cells[jRow * mesh->Nx].cy);

    ts << makeHeader(params, title, snap.iter, snap.simTime);
    ts << "i,x,rho,u,v,p,mach\n";

    const int Nx = mesh->Nx;

    for(int i = 0; i < Nx; ++i)
    {
        int idx = jRow * Nx + i;

        const Cell&       cell = mesh->cells[idx];
        const EulerState& U    = snap.cells[idx];

        double rho = U.rho;
        double u   = (rho > 1e-30) ? U.rho_u / rho : 0.0;
        double v   = (rho > 1e-30) ? U.rho_v / rho : 0.0;
        double p   = EulerPhysics::pressure(U);
        double a   = EulerPhysics::soundSpeed(U);
        double spd = std::sqrt(u*u + v*v);
        double M   = (a > 1e-30) ? spd / a : 0.0;

        ts << i       << ","
           << cell.cx << ","
           << rho      << ","
           << u        << ","
           << v        << ","
           << p        << ","
           << M        << "\n";
    }

    return true;
}
