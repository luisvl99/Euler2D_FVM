#pragma once

#include <QString>
#include <vector>

#include "io_csv.h"
#include "mesh_structured.h"
#include "physics_eulerstate.h"

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
//  Qt side of the CSV export.  The file format lives in io_csv.h (shared
//  with the command-line runner); this class opens the file (Unicode paths
//  included) and turns failures into messages for a QMessageBox.
//
//  All functions return true on success.  On failure they return false and,
//  if errorOut != nullptr, set it to a human-readable message.
// ---------------------------------------------------------------------------
class ExportManager
{
public:

    static bool writeResiduals(const QString&             path,
                               const RunInfo&             info,
                               const std::vector<double>& residuals,
                               const std::vector<double>& dtHistory,
                               QString*                   errorOut = nullptr);

    static bool writeSnapshot(const QString&        path,
                              const RunInfo&        info,
                              const StructuredMesh& mesh,
                              const CellSnapshot&   snap,
                              QString*              errorOut = nullptr);

    static bool writeLineProbe(const QString&        path,
                               const RunInfo&        info,
                               const StructuredMesh& mesh,
                               const CellSnapshot&   snap,
                               int                   jRow,
                               QString*              errorOut = nullptr);
};
