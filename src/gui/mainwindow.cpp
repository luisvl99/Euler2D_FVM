#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "mesh_structured.h"
#include "physics_eulerphysics.h"
#include "gui_exportmanager.h"

#include <QAction>
#include <QFileDialog>
#include <QLabel>
#include <QMessageBox>
#include <QSpinBox>
#include <QLineEdit>
#include <QPushButton>

#include <algorithm>
#include <cmath>
#include <exception>
#include <iostream>
#include <QString>

// ===========================================================================
//  Constructor / Destructor
// ===========================================================================

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    setWindowState(Qt::WindowMaximized);

    // -----------------------------------------------------------------------
    // Checkbox → MeshViewer connections (unchanged)
    // -----------------------------------------------------------------------
    connect(ui->checkDrawCenters,
            &QCheckBox::toggled,
            ui->meshViewer,
            &MeshViewer::setDrawCellCenters);

    connect(ui->checkDrawNormals,
            &QCheckBox::toggled,
            ui->meshViewer,
            &MeshViewer::setDrawFaceNormals);

    // -----------------------------------------------------------------------
    // Field combo box
    // -----------------------------------------------------------------------
    ui->comboField->addItem("None (wireframe)",   (int)FieldType::None);
    ui->comboField->addItem("Density ρ",          (int)FieldType::Density);
    ui->comboField->addItem("Velocity |V|",       (int)FieldType::VelocityMag);
    ui->comboField->addItem("Velocity u",         (int)FieldType::VelocityU);
    ui->comboField->addItem("Velocity v",         (int)FieldType::VelocityV);
    ui->comboField->addItem("Pressure p",         (int)FieldType::Pressure);
    ui->comboField->addItem("Mach number",        (int)FieldType::MachNumber);
    ui->comboField->addItem("Total energy E",     (int)FieldType::TotalEnergy);
    ui->comboField->addItem("Residual |R|",       (int)FieldType::Residual);

    // -----------------------------------------------------------------------
    // Scheme combo box
    // -----------------------------------------------------------------------
    ui->comboScheme->addItem("SSP-RK2",        QString("RK2"));
    ui->comboScheme->addItem("Forward Euler",  QString("ForwardEuler"));

    ui->comboFluxScheme->addItem("Rusanov",  QString("Rusanov"));
    ui->comboFluxScheme->addItem("HLLC",     QString("HLLC"));

    ui->comboReconstruction->addItem("Piecewise constant", QString("PC"));
    ui->comboReconstruction->addItem("MUSCL (MinMod)",     QString("MUSCL"));

    ui->comboCase->addItem("Channel flow",           (int)TestCase::Channel);
    ui->comboCase->addItem("Sod shock tube",         (int)TestCase::Sod);
    ui->comboCase->addItem("Shu-Osher",              (int)TestCase::ShuOsher);
    ui->comboCase->addItem("2-D Riemann, config 3",  (int)TestCase::Riemann3);
    ui->comboCase->addItem("2-D Riemann, config 4",  (int)TestCase::Riemann4);
    ui->comboCase->addItem("2-D Riemann, config 6",  (int)TestCase::Riemann6);
    ui->comboCase->addItem("2-D Riemann, config 12", (int)TestCase::Riemann12);

    // -----------------------------------------------------------------------
    // Export menu
    // -----------------------------------------------------------------------
    connect(ui->checkExportProbe, &QCheckBox::toggled,
            ui->spinExportJ,      &QSpinBox::setEnabled);


    // -----------------------------------------------------------------------
    // Solver timer
    //
    //  Interval = 0 ms means the timer fires as soon as the event loop is
    //  idle — i.e. as fast as possible while still processing GUI events
    //  (mouse clicks, repaints, etc.) between ticks.
    // -----------------------------------------------------------------------
    m_solverTimer = new QTimer(this);
    m_solverTimer->setInterval(0);
    connect(m_solverTimer, &QTimer::timeout,
            this,           &MainWindow::onSolverStep);

    //Residual plot
    m_residualPlotter = new ResidualPlotter(this);
    m_residualPlotter->setWindowFlags(Qt::Window);
    m_residualPlotter->setWindowTitle("Residual History");
    m_residualPlotter->resize(600, 350);

    // Stop button starts disabled (nothing running yet)
    ui->btnStop->setEnabled(false);
    ui->lblStatus->setText("Ready.");
}

MainWindow::~MainWindow()
{
    // Stop timer before destroying solver / mesh
    if(m_solverTimer->isActive())
        m_solverTimer->stop();

    delete ui;
}

// ===========================================================================
//  Mesh generation
// ===========================================================================

void MainWindow::on_btnGenerateMesh_clicked()
{
    // If solver is running, stop it first
    if(m_solverTimer->isActive())
        m_solverTimer->stop();

    int Nx = ui->spinNx->value();
    int Ny = ui->spinNy->value();

    double Lx = ui->spinLx->value();
    double Ly = ui->spinLy->value();

    mesh   = std::make_unique<StructuredMesh>(Nx, Ny, Lx, Ly);
    solver = nullptr;         // old solver is invalid once mesh changes
    m_snapshots.clear();      // old snapshots no longer match the cell count
    ui->spinExportJ->setMaximum(Ny - 1);

    std::cout << "[Mesh] Generated " << Nx << "×" << Ny
              << "  (" << mesh->cells.size() << " cells)\n";

    ui->meshViewer->setMesh(mesh.get());
    setControlsEnabled(false);   // not running
    ui->lblStatus->setText(
        QString("Mesh ready: %1×%2 (%3 cells)")
            .arg(Nx).arg(Ny).arg(mesh->cells.size()));
}

// ===========================================================================
//  Run button — initialises solver and starts the timer loop
// ===========================================================================

void MainWindow::on_btnRunSolver_clicked()
{
    if(!mesh)
    {
        std::cerr << "[MainWindow] Error: generate a mesh first.\n";
        ui->lblStatus->setText("Generate a mesh first.");
        return;
    }

    const TestCase tc = static_cast<TestCase>(ui->comboCase->currentData().toInt());

    // Selecting Shu-Osher only pre-fills the spin boxes; if the mesh was not
    // regenerated, its setup lands on the wrong domain
    if(tc == TestCase::ShuOsher && std::abs(mesh->Lx - SHU_OSHER_LX) > 1e-9)
    {
        const auto answer = QMessageBox::question(
            this, "Shu-Osher domain",
            QString("Shu-Osher is defined on x in [-5, 5], so it needs Lx = %1.\n"
                    "The current mesh has Lx = %2. Click Generate Mesh after "
                    "selecting the case to use the canonical mesh.\n\n"
                    "Run anyway?")
                .arg(SHU_OSHER_LX).arg(mesh->Lx),
            QMessageBox::Yes | QMessageBox::No, QMessageBox::No);

        if(answer != QMessageBox::Yes)
        {
            ui->lblStatus->setText(
                QString("Run cancelled: Shu-Osher needs a mesh with Lx = %1.").arg(SHU_OSHER_LX));
            return;
        }
    }

    // ------------------------------------------------------------------
    // Flow parameters (Channel uses them; the other cases set their own)
    // ------------------------------------------------------------------
    FlowParameters fp;
    fp.rho_inf = ui->spinRho->value();
    fp.u_inf   = ui->spinU->value();
    fp.v_inf   = ui->spinV->value();
    fp.p_inf   = ui->spinP->value();
    fp.p_back  = ui->spinPback->value();

    // ------------------------------------------------------------------
    // Boundary tags + initial condition — must happen BEFORE building the
    // solver (its constructor registers one BC object per BoundaryType)
    // ------------------------------------------------------------------
    const double tEnd = applyTestCase(tc, *mesh, fp, ui->checkSubsonic->isChecked());

    m_simParams = buildSimParams(fp);

    // ------------------------------------------------------------------
    // Build solver
    // ------------------------------------------------------------------
    solver = std::make_unique<EulerSolver>(mesh.get(), fp);
    solver->setCFL(m_simParams.CFL);
    solver->setTimeScheme(
        m_simParams.scheme == "RK2"
            ? EulerSolver::TimeScheme::RK2
            : EulerSolver::TimeScheme::ForwardEuler);

    solver->setFluxScheme(
        ui->comboFluxScheme->currentData().toString() == "HLLC"
            ? EulerSolver::FluxScheme::HLLC
            : EulerSolver::FluxScheme::Rusanov);

    solver->setReconstructionScheme(
        ui->comboReconstruction->currentData().toString() == "MUSCL"
            ? EulerSolver::ReconstructionScheme::MUSCL
            : EulerSolver::ReconstructionScheme::PiecewiseConstant);

    solver->setTEnd(tEnd);

    // ------------------------------------------------------------------
    // Reset counters and snapshots
    // ------------------------------------------------------------------
    m_snapshots.clear();
    m_snapshotEvery  = ui->spinSnapshotEvery->value();
    m_maxIter        = ui->spinMaxIter->value();
    m_stepsPerFrame  = ui->spinStepsPerFrame->value();
    m_currentIter    = 0;
    storeSnapshot(0, 0.0);

    // ------------------------------------------------------------------
    // 1-D cases: live profile viewer (exact overlay for Sod)
    // ------------------------------------------------------------------
    delete m_sodViewer;
    m_sodViewer = nullptr;

    if(tc == TestCase::Sod || tc == TestCase::ShuOsher)
    {
        m_sodViewer = new SodViewer(this);
        m_sodViewer->setWindowFlags(Qt::Window);
        m_sodViewer->resize(740, 660);
        // Middle j-row is the probe line (quasi-1D anyway)
        m_sodViewer->configure(ui->comboCase->currentText(), mesh->Lx, mesh->Ny / 2);
        if(tc == TestCase::Sod)
            m_sodViewer->setExact({1.0, 0.0, 1.0}, {0.125, 0.0, 0.1}, mesh->Lx * 0.5);
        m_sodViewer->show();
        m_sodViewer->raise();
    }

    // ------------------------------------------------------------------
    // Start
    // ------------------------------------------------------------------
    setControlsEnabled(true);
    ui->lblStatus->setText(QString("Running (%1)...").arg(ui->comboCase->currentText()));

    m_residualPlotter->setData({});
    m_residualPlotter->show();
    m_residualPlotter->raise();
    m_solverTimer->start();
}

// ===========================================================================
//  Stop button
// ===========================================================================

void MainWindow::on_btnStop_clicked()
{
    if(m_solverTimer->isActive())
    {
        m_solverTimer->stop();
        setControlsEnabled(false);

        // Show final status
        if(solver)
        {
            std::cout << "[MainWindow] Solver stopped by user at iter "
                      << m_currentIter << "\n";
        }
        ui->lblStatus->setText(
            QString("Stopped at iter %1 / %2")
                .arg(m_currentIter).arg(m_maxIter));
    }
}

// ===========================================================================
//  Timer slot — called every event-loop cycle while the solver is running
//
//  Advances the solution by m_stepsPerFrame steps, then updates the viewer.
//  Stops automatically when the iteration count reaches m_maxIter.
// ===========================================================================

void MainWindow::onSolverStep()
{
    if(!solver) { m_solverTimer->stop(); return; }

    // step() throws when the solution becomes non-physical (blow-up)
    QString failure;
    try
    {
        for(int s = 0; s < m_stepsPerFrame
                      && m_currentIter < m_maxIter
                      && solver->simTime() < solver->tEnd(); ++s)
        {
            solver->step();
            ++m_currentIter;

            // Store snapshot at the configured interval
            if(m_snapshotEvery > 0 && m_currentIter % m_snapshotEvery == 0)
                storeSnapshot(m_currentIter, solver->simTime());
        }
    }
    catch(const std::exception& e)
    {
        failure = QString::fromUtf8(e.what());
    }

    ui->meshViewer->update();
    m_residualPlotter->setData(solver->residualHistory());
    if(m_sodViewer)
        m_sodViewer->refresh(mesh.get(), solver->simTime());
    updateStatusLabel();

    if(!failure.isEmpty())
    {
        // The failed state stays on screen so the user can see where it started
        m_solverTimer->stop();
        setControlsEnabled(false);
        ui->lblStatus->setText(
            QString("Stopped: non-physical state at iter %1").arg(m_currentIter + 1));
        QMessageBox::critical(this, "Solver stopped", failure);
        return;
    }

    if(m_currentIter >= m_maxIter || solver->simTime() >= solver->tEnd())
    {
        // Always capture a final snapshot so the last state is exportable
        if(m_snapshots.empty() ||
            m_snapshots.back().iter != m_currentIter)
        {
            storeSnapshot(m_currentIter, solver->simTime());
        }

        m_solverTimer->stop();
        setControlsEnabled(false);
        ui->lblStatus->setText(
            QString("Finished: %1 iters  t=%2  ||R||=%3")
                .arg(m_currentIter)
                .arg(solver->simTime(),               0, 'e', 4)
                .arg(solver->residualHistory().back(), 0, 'e', 3));

        std::cout << "[MainWindow] Finished after " << m_currentIter << " iterations.\n";
    }
}









// ===========================================================================
//  Field combo box
// ===========================================================================

void MainWindow::on_comboField_currentIndexChanged(int index)
{
    FieldType ft = static_cast<FieldType>(
        ui->comboField->itemData(index).toInt());

    ui->meshViewer->setActiveField(ft);
}

// ===========================================================================
//  Case combo box — pre-fill the canonical mesh of the selected case
//  (the user still has to press Generate Mesh, and may change the values)
// ===========================================================================

void MainWindow::on_comboCase_currentIndexChanged(int index)
{
    const TestCase tc = static_cast<TestCase>(ui->comboCase->itemData(index).toInt());
    if(tc == TestCase::Channel || tc == TestCase::Sod) return;   // any domain works

    const bool shuOsher = (tc == TestCase::ShuOsher);
    ui->spinNx->setValue(shuOsher ? 400  : 200);
    ui->spinNy->setValue(shuOsher ? 1    : 200);
    ui->spinLx->setValue(shuOsher ? SHU_OSHER_LX : 1.0);
    ui->spinLy->setValue(1.0);

    // Let the run stop at the case's end time, not at the iteration cap
    ui->spinMaxIter->setValue(std::max(ui->spinMaxIter->value(), 100000));
    ui->comboField->setCurrentIndex(ui->comboField->findData((int)FieldType::Density));
}

void MainWindow::on_btnExportBrowse_clicked()
{
    QString d = QFileDialog::getExistingDirectory(
        this, "Choose export directory", ui->lineExportDir->text());
    if(!d.isEmpty())
        ui->lineExportDir->setText(d);
}

void MainWindow::on_btnExport_clicked()
{
    QString dir = ui->lineExportDir->text().trimmed();
    if(dir.isEmpty())
    {
        QMessageBox::warning(this, "Export", "Please choose an output directory first.");
        return;
    }

    int written = 0, failed = 0;
    QStringList errors;

    if(ui->checkExportResiduals->isChecked())
    {
        if(!solver || solver->residualHistory().empty())
        {
            errors << "No residual data (run the solver first).";
            ++failed;
        }
        else
        {
            QString err;
            QString path = dir + "/residuals.csv";
            if(ExportManager::writeResiduals(path, solver->residualHistory(),
                                              solver->dtHistory(), m_simParams, &err))
                ++written;
            else { ++failed; errors << err; }
        }
    }

    if(ui->checkExportSnapshots->isChecked())
    {
        if(m_snapshots.empty())
        {
            errors << "No snapshots available.";
            ++failed;
        }
        else
        {
            for(const CellSnapshot& snap : m_snapshots)
            {
                QString path = QString("%1/snapshot_iter_%2.csv")
                .arg(dir)
                    .arg(snap.iter, 7, 10, QChar('0'));
                QString err;
                if(ExportManager::writeSnapshot(path, mesh.get(), snap, m_simParams, &err))
                    ++written;
                else { ++failed; errors << err; }
            }
        }
    }

    if(ui->checkExportProbe->isChecked())
    {
        if(m_snapshots.empty())
        {
            errors << "No snapshots available.";
            ++failed;
        }
        else
        {
            int jRow = ui->spinExportJ->value();
            for(const CellSnapshot& snap : m_snapshots)
            {
                QString path = QString("%1/probe_j%2_iter_%3.csv")
                .arg(dir)
                    .arg(jRow,      2, 10, QChar('0'))
                    .arg(snap.iter, 7, 10, QChar('0'));
                QString err;
                if(ExportManager::writeLineProbe(path, mesh.get(), snap, jRow,
                                                  m_simParams, &err))
                    ++written;
                else { ++failed; errors << err; }
            }
        }
    }

    if(failed == 0)
        QMessageBox::information(this, "Export",
                                 QString("Done. %1 file(s) written to:\n%2").arg(written).arg(dir));
    else
        QMessageBox::warning(this, "Export",
                             QString("%1 file(s) written, %2 failed.\n\n%3")
                                 .arg(written).arg(failed).arg(errors.join("\n")));
}


// ===========================================================================
//  Helpers
// ===========================================================================

void MainWindow::storeSnapshot(int iter, double simTime)
{
    if(!mesh) return;

    CellSnapshot snap;
    snap.iter    = iter;
    snap.simTime = simTime;
    snap.cells.resize(mesh->cells.size());

    for(std::size_t k = 0; k < mesh->cells.size(); ++k)
        snap.cells[k] = mesh->cells[k].U;

    m_snapshots.push_back(std::move(snap));
}

// Export metadata for the run being started.  Mesh and flow data come from
// the objects actually used, not from the spin boxes: picking a case changes
// the mesh spin boxes without regenerating the mesh, and every case except
// Channel replaces the GUI flow state with its own.
SimParameters MainWindow::buildSimParams(const FlowParameters& fp) const
{
    SimParameters p;
    p.Nx  = mesh->Nx;
    p.Ny  = mesh->Ny;
    p.Lx  = mesh->Lx;
    p.Ly  = mesh->Ly;

    p.rho_inf = fp.rho_inf;
    p.u_inf   = fp.u_inf;
    p.v_inf   = fp.v_inf;
    p.p_inf   = fp.p_inf;

    p.CFL           = ui->spinCFL->value();
    p.scheme        = ui->comboScheme->currentData().toString();
    p.maxIter       = ui->spinMaxIter->value();
    p.snapshotEvery = ui->spinSnapshotEvery->value();

    p.fluxScheme    = ui->comboFluxScheme->currentData().toString();
    p.reconstruction = ui->comboReconstruction->currentData().toString();

    return p;
}

void MainWindow::setControlsEnabled(bool running)
{
    ui->btnGenerateMesh->setEnabled(!running);
    ui->btnRunSolver->setEnabled(!running);
    ui->spinNx->setEnabled(!running);
    ui->spinNy->setEnabled(!running);
    ui->spinLx->setEnabled(!running);
    ui->spinLy->setEnabled(!running);
    ui->spinRho->setEnabled(!running);
    ui->spinU->setEnabled(!running);
    ui->spinV->setEnabled(!running);
    ui->spinP->setEnabled(!running);
    ui->spinCFL->setEnabled(!running);
    ui->spinMaxIter->setEnabled(!running);
    ui->spinStepsPerFrame->setEnabled(!running);
    ui->spinSnapshotEvery->setEnabled(!running);
    ui->comboScheme->setEnabled(!running);
    ui->comboFluxScheme->setEnabled(!running);
    ui->comboReconstruction->setEnabled(!running);
    ui->checkSubsonic->setEnabled(!running);
    ui->comboCase->setEnabled(!running);
    ui->spinPback->setEnabled(!running);
    ui->btnExport->setEnabled(!running);
    ui->btnStop->setEnabled(running);
}

void MainWindow::updateStatusLabel()
{
    if(!solver || solver->residualHistory().empty()) return;

    ui->lblStatus->setText(
        QString("Iter %1/%2  t=%3  dt=%4  ||R||=%5  snaps=%6")
            .arg(m_currentIter)
            .arg(m_maxIter)
            .arg(solver->simTime(),               0, 'e', 3)
            .arg(solver->lastDt(),                0, 'e', 3)
            .arg(solver->residualHistory().back(), 0, 'e', 3)
            .arg(m_snapshots.size()));
}

