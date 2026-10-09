#include "gui_exportmanager.h"

#include <filesystem>
#include <fstream>

// ===========================================================================
//  writeFile — opens path and lets write(out) fill it.
//  Returns false (sets errorOut) if the file cannot be opened or written.
// ===========================================================================
template<class Write>
static bool writeFile(const QString& path, QString* errorOut, Write write)
{
    std::ofstream out(std::filesystem::path(path.toStdWString()));
    if(out)
        write(out);
    out.close();

    if(!out)
    {
        if(errorOut)
            *errorOut = QString("Cannot write file:\n%1").arg(path);
        return false;
    }
    return true;
}

bool ExportManager::writeResiduals(const QString&             path,
                                   const RunInfo&             info,
                                   const std::vector<double>& residuals,
                                   const std::vector<double>& dtHistory,
                                   QString*                   errorOut)
{
    if(residuals.empty())
    {
        if(errorOut) *errorOut = "No residual data to export (run the solver first).";
        return false;
    }

    return writeFile(path, errorOut, [&](std::ostream& out) {
        writeResidualsCsv(out, info, residuals, dtHistory);
    });
}

bool ExportManager::writeSnapshot(const QString&        path,
                                  const RunInfo&        info,
                                  const StructuredMesh& mesh,
                                  const CellSnapshot&   snap,
                                  QString*              errorOut)
{
    if(snap.cells.size() != mesh.cells.size())
    {
        if(errorOut) *errorOut = "No snapshot data available for this mesh.";
        return false;
    }

    return writeFile(path, errorOut, [&](std::ostream& out) {
        writeSnapshotCsv(out, info, mesh, snap.cells, snap.iter, snap.simTime);
    });
}

bool ExportManager::writeLineProbe(const QString&        path,
                                   const RunInfo&        info,
                                   const StructuredMesh& mesh,
                                   const CellSnapshot&   snap,
                                   int                   jRow,
                                   QString*              errorOut)
{
    if(snap.cells.size() != mesh.cells.size())
    {
        if(errorOut) *errorOut = "No snapshot data available for this mesh.";
        return false;
    }
    if(jRow < 0 || jRow >= mesh.Ny)
    {
        if(errorOut) *errorOut = QString("j=%1 is out of range [0, %2].")
                                     .arg(jRow).arg(mesh.Ny - 1);
        return false;
    }

    return writeFile(path, errorOut, [&](std::ostream& out) {
        writeLineProbeCsv(out, info, mesh, snap.cells, jRow, snap.iter, snap.simTime);
    });
}
