#pragma once
#include <QWidget>
#include <vector>
#include "sod_exact.h"

class StructuredMesh;

// ---------------------------------------------------------------------------
// SodViewer
//
//  Floating window that shows 3 stacked panels (ρ, u, p vs x) updated live
//  during a 1-D run (Sod, Shu-Osher).
//
//  Each panel draws:
//    - Red continuous line  : exact Riemann solution (only after setExact())
//    - Blue filled circles  : FVM cell-centre values at j = m_jRow
//
//  Typical usage in MainWindow
//  ---------------------------
//    // Once, when the run is launched:
//    m_sodViewer->configure("Sod shock tube", mesh->Lx, Ny/2);
//    m_sodViewer->setExact(leftState, rightState, x0);   // optional
//    m_sodViewer->show();
//
//    // Every timer tick (inside onSolverStep):
//    m_sodViewer->refresh(mesh.get(), solver->simTime());
// ---------------------------------------------------------------------------

class SodViewer : public QWidget
{
    Q_OBJECT

public:
    explicit SodViewer(QWidget* parent = nullptr);

    // Call once per run: window title, domain length and j-row to probe
    void configure(const QString& title, double Lx, int jRow);

    // Optional: overlay the exact Riemann solution for states L | R at x0
    void setExact(SodExact::State L, SodExact::State R, double x0);

    // Call every solver step — extracts FVM data and recomputes exact line
    void refresh(const StructuredMesh* mesh, double t);

protected:
    void paintEvent(QPaintEvent*) override;

private:
    // ---- Exact solver ---------------------------------------------------
    SodExact m_exact;
    bool     m_hasExact   = false;
    bool     m_configured = false;
    QString  m_title;

    double m_x0  = 0.5;
    double m_Lx  = 1.0;
    int    m_jRow = 0;
    double m_time = 0.0;

    // ---- Cached data for painting (filled by refresh()) -----------------
    struct FvmPt { double x, rho, u, p; };
    std::vector<FvmPt>        m_fvm;
    std::vector<SodExact::Point> m_exactPts;   // 500-point exact line

    // ---- Layout helpers -------------------------------------------------
    // Returns the inner plot rectangle for panel index (0=rho, 1=u, 2=p)
    // given the current widget dimensions.
    struct PlotRect { int left, top, w, h; };
    PlotRect plotRect(int panelIdx) const;

    // Draw one panel — varIdx selects the variable
    void drawPanel(QPainter& painter,
                   const PlotRect& pr,
                   int varIdx,
                   const QString& varLabel) const;
};
