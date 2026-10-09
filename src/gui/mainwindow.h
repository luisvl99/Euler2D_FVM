#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QTimer>
#include <memory>

#include "mesh_structured.h"
#include "gui_meshviewer.h"
#include "gui_fieldtype.h"
#include "physics_solver.h"
#include "gui_exportmanager.h"
#include "gui_residualplotter.h"
#include "physics_testcases.h"
#include "gui_sodviewer.h"


QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private:
    // -----------------------------------------------------------------------
    // UI
    // -----------------------------------------------------------------------
    Ui::MainWindow *ui;

    // -----------------------------------------------------------------------
    // Simulation objects
    // -----------------------------------------------------------------------
    std::unique_ptr<StructuredMesh> mesh;
    std::unique_ptr<EulerSolver>    solver;

    // -----------------------------------------------------------------------
    // Timer-driven solver loop
    //
    //  m_solverTimer fires at 0 ms (as fast as the event loop allows).
    //  Each tick advances the solution by m_stepsPerFrame steps and then
    //  requests a repaint.  The loop stops automatically when m_currentIter
    //  reaches m_maxIter, or immediately when the user presses Stop.
    // -----------------------------------------------------------------------
    QTimer* m_solverTimer  = nullptr;
    int     m_currentIter  = 0;   // steps completed since last Run
    int     m_maxIter      = 0;   // target set by spinMaxIter at Run time
    int     m_stepsPerFrame = 1;  // steps per timer tick (from spinStepsPerFrame)


    // -----------------------------------------------------------------------
    // Snapshot storage
    //
    //  Filled during the run every spinSnapshotEvery iterations.
    //  Passed to ExportManager on demand — the live mesh is NOT used
    //  by the exporter so there is no risk of exporting mid-step state.
    // -----------------------------------------------------------------------
    std::vector<CellSnapshot> m_snapshots;
    int                       m_snapshotEvery = 100;

    // Export metadata of the current run (set when Run is clicked)
    RunInfo m_runInfo;

    //Residual plotter
    ResidualPlotter* m_residualPlotter = nullptr;
    SodViewer*       m_sodViewer       = nullptr;   // 1-D cases only

    // -----------------------------------------------------------------------
    // Helpers
    // -----------------------------------------------------------------------
    void setControlsEnabled(bool running);  // grey-out UI while solver runs
    void updateStatusLabel();               // refresh lblStatus text

    // Capture a snapshot of the current cell states
    void storeSnapshot(int iter, double simTime);



private slots:
    // Qt Designer auto-connect (name matches widget objectName)
    void on_btnGenerateMesh_clicked();
    void on_btnRunSolver_clicked();
    void on_btnStop_clicked();
    void on_comboField_currentIndexChanged(int index);
    void on_comboCase_currentIndexChanged(int index);

    void on_btnExportBrowse_clicked();
    void on_btnExport_clicked();

    // Timer callback — NOT auto-connected, connected manually in constructor
    void onSolverStep();
};

#endif // MAINWINDOW_H
