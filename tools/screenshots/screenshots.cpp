// ---------------------------------------------------------------------------
// GUI screenshots for the README
//
//  Drives the real MainWindow through its widgets (as a user would), waits
//  for each run to finish and saves window grabs to the folder given as the
//  first argument.  See CMakeLists.txt next to this file for how to build
//  and run it.
// ---------------------------------------------------------------------------

#include <QApplication>
#include <QCheckBox>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QElapsedTimer>
#include <QFont>
#include <QPushButton>
#include <QSpinBox>

#include "gui_residualplotter.h"
#include "gui_sodviewer.h"
#include "mainwindow.h"

template<class T>
static T* child(QWidget& w, const char* name)
{
    T* c = w.findChild<T*>(name);
    if(!c) qFatal("missing widget %s", name);
    return c;
}

static void select(QWidget& w, const char* combo, const QVariant& data)
{
    QComboBox* box = child<QComboBox>(w, combo);
    const int i = box->findData(data);
    if(i < 0) qFatal("no item %s in %s", qPrintable(data.toString()), combo);
    box->setCurrentIndex(i);
}

static void mesh(QWidget& w, int nx, int ny, double lx, double ly)
{
    child<QSpinBox>(w, "spinNx")->setValue(nx);
    child<QSpinBox>(w, "spinNy")->setValue(ny);
    child<QDoubleSpinBox>(w, "spinLx")->setValue(lx);
    child<QDoubleSpinBox>(w, "spinLy")->setValue(ly);
    child<QPushButton>(w, "btnGenerateMesh")->click();
}

// Clicks Run and processes events until the run has finished
static void runToEnd(QApplication& app, QWidget& w)
{
    QPushButton* run = child<QPushButton>(w, "btnRunSolver");
    run->click();
    QElapsedTimer t;
    t.start();
    while(!run->isEnabled() && t.elapsed() < 600000)
        app.processEvents(QEventLoop::AllEvents, 50);
    for(int k = 0; k < 10; ++k)
        app.processEvents();
}

int main(int argc, char** argv)
{
    QApplication app(argc, argv);
    QApplication::setFont(QFont("Segoe UI", 9));
    const QString out = argc > 1 ? argv[1] : ".";

    MainWindow w;
    w.setWindowState(Qt::WindowNoState);
    w.resize(1500, 880);
    w.show();
    app.processEvents();

    select(w, "comboFluxScheme", QString("HLLC"));
    select(w, "comboReconstruction", QString("MUSCL"));
    select(w, "comboScheme", QString("RK2"));

    // 1. Main window: 2-D Riemann problem, configuration 3, density
    select(w, "comboCase", int(TestCase::Riemann3));
    child<QSpinBox>(w, "spinStepsPerFrame")->setValue(50);
    mesh(w, 200, 200, 1.0, 1.0);
    runToEnd(app, w);
    w.grab().save(out + "/gui_riemann.png");

    // 2. Sod shock tube: live profile viewer with the exact solution
    select(w, "comboCase", int(TestCase::Sod));
    child<QSpinBox>(w, "spinStepsPerFrame")->setValue(5);
    mesh(w, 200, 1, 1.0, 0.1);
    runToEnd(app, w);
    if(auto* s = w.findChild<SodViewer*>())
        s->grab().save(out + "/gui_sod_profiles.png");

    // 3. Subsonic channel flow from rest: residual history to steady state
    select(w, "comboCase", int(TestCase::Channel));
    child<QDoubleSpinBox>(w, "spinU")->setValue(68.0);
    child<QCheckBox>(w, "checkSubsonic")->setChecked(true);
    child<QSpinBox>(w, "spinMaxIter")->setValue(10000);
    child<QSpinBox>(w, "spinStepsPerFrame")->setValue(100);
    select(w, "comboField", int(FieldType::MachNumber));
    mesh(w, 100, 20, 1.0, 0.2);
    runToEnd(app, w);
    if(auto* r = w.findChild<ResidualPlotter*>())
        r->grab().save(out + "/gui_residuals.png");

    return 0;
}
