#include "gui_sodviewer.h"
#include "mesh_structured.h"
#include "physics_eulerphysics.h"

#include <QPainter>
#include <QPen>
#include <QFont>
#include <QFontMetrics>

#include <algorithm>
#include <cmath>
#include <limits>

// ===========================================================================
//  Constructor
// ===========================================================================

SodViewer::SodViewer(QWidget* parent)
    : QWidget(parent)
{
    setMinimumSize(720, 620);
}

// ===========================================================================
//  configure  —  call once per run before showing the window
// ===========================================================================

void SodViewer::configure(const QString& title, double Lx, int jRow)
{
    m_title = title;
    setWindowTitle(title);
    m_Lx   = Lx;
    m_jRow = jRow;
    m_time = 0.0;
    m_fvm.clear();
    m_exactPts.clear();
    m_hasExact   = false;
    m_configured = true;
    update();
}

void SodViewer::setExact(SodExact::State L, SodExact::State R, double x0)
{
    m_exact.setStates(L, R, 1.4);
    m_x0       = x0;
    m_hasExact = true;
}

// ===========================================================================
//  refresh  —  called every solver step
// ===========================================================================

void SodViewer::refresh(const StructuredMesh* mesh, double t)
{
    if(!mesh || !m_configured) return;

    m_time = t;

    // ---- Extract FVM primitives at j = m_jRow ---------------------------
    m_fvm.clear();
    m_fvm.reserve(mesh->Nx);

    for(int i = 0; i < mesh->Nx; ++i)
    {
        const int   idx  = m_jRow * mesh->Nx + i;
        const auto& cell = mesh->cells[idx];
        const auto& U    = cell.U;

        if(U.rho < 1e-20) continue;    // skip degenerate cells

        FvmPt pt;
        pt.x   = cell.cx;
        pt.rho = U.rho;
        pt.u   = U.rho_u / U.rho;
        pt.p   = EulerPhysics::pressure(U);
        m_fvm.push_back(pt);
    }

    // ---- Recompute exact solution at 500 points -------------------------
    if(m_hasExact)
        m_exactPts = m_exact.sampleLine(500, m_Lx, t, m_x0);

    update();   // trigger repaint
}

// ===========================================================================
//  plotRect  —  inner plot area for panel panelIdx (0, 1, 2)
//
//  Layout (pixels):
//    Title strip        : 32 px at top
//    Left margin        : 60 px  (y-axis numbers)
//    Right margin       : 20 px
//    Bottom strip       : 34 px  (x-axis label, shown on lowest panel)
//    Between panels     : 8 px gap
// ===========================================================================

SodViewer::PlotRect SodViewer::plotRect(int panelIdx) const
{
    const int titleH  = 32;
    const int bottomH = 34;
    const int ml      = 60;
    const int mr      = 20;
    const int gap     = 8;

    const int totalPlotH = height() - titleH - bottomH - 2 * gap;
    const int panelH     = totalPlotH / 3;
    const int panelW     = width() - ml - mr;

    const int top = titleH + panelIdx * (panelH + gap);

    return { ml, top, panelW, panelH };
}

// ===========================================================================
//  paintEvent
// ===========================================================================

void SodViewer::paintEvent(QPaintEvent*)
{
    QPainter painter(this);
    painter.setRenderHint(QPainter::Antialiasing);

    // Background
    painter.fillRect(rect(), QColor(24, 24, 28));

    // Title
    {
        QFont f = painter.font();
        f.setBold(true);
        f.setPointSize(10);
        painter.setFont(f);
        painter.setPen(QColor(220, 220, 220));

        QString title = m_configured
                            ? QString("%1  ·  t = %2").arg(m_title).arg(m_time, 0, 'f', 4)
                            : QString("(not configured)");
        painter.drawText(QRect(0, 4, width(), 26),
                         Qt::AlignHCenter | Qt::AlignVCenter, title);
    }

    if(!m_configured)
    {
        painter.setPen(QColor(160, 160, 160));
        painter.drawText(rect(), Qt::AlignCenter,
                         "Call configure() before the first run.");
        return;
    }

    // Draw the three panels
    drawPanel(painter, plotRect(0), 0, "rho  (density)");
    drawPanel(painter, plotRect(1), 1, "u  (velocity)");
    drawPanel(painter, plotRect(2), 2, "p  (pressure)");

    // x-axis label at the very bottom
    {
        QFont f = painter.font();
        f.setBold(false);
        f.setPointSize(8);
        painter.setFont(f);
        painter.setPen(QColor(160, 160, 160));
        painter.drawText(QRect(plotRect(2).left, height() - 22,
                               plotRect(2).w, 20),
                         Qt::AlignHCenter, "x");
    }
}

// ===========================================================================
//  drawPanel
//
//  varIdx:  0 = rho,  1 = u,  2 = p
// ===========================================================================

void SodViewer::drawPanel(QPainter& painter, const PlotRect& pr,
                          int varIdx, const QString& varLabel) const
{
    // ---- Lambdas to extract variable from cached data -------------------
    auto getExact = [varIdx](const SodExact::Point& pt) -> double {
        switch(varIdx) {
        case 0: return pt.rho;
        case 1: return pt.u;
        default: return pt.p;
        }
    };
    auto getFvm = [varIdx](const FvmPt& pt) -> double {
        switch(varIdx) {
        case 0: return pt.rho;
        case 1: return pt.u;
        default: return pt.p;
        }
    };

    // ---- Compute y range (from combined exact + FVM) --------------------
    double yMin =  std::numeric_limits<double>::max();
    double yMax = -std::numeric_limits<double>::max();

    for(const auto& pt : m_exactPts)
    {
        double v = getExact(pt);
        yMin = std::min(yMin, v);
        yMax = std::max(yMax, v);
    }
    for(const auto& pt : m_fvm)
    {
        double v = getFvm(pt);
        yMin = std::min(yMin, v);
        yMax = std::max(yMax, v);
    }

    // Fallback when no data
    if(yMin > yMax) { yMin = 0.0; yMax = 1.0; }

    // Ensure a non-degenerate range with 8% padding
    const double span = (yMax - yMin > 1e-12) ? (yMax - yMin) : 0.2;
    yMin -= 0.08 * span;
    yMax += 0.08 * span;

    const double xMin = 0.0;
    const double xMax = m_Lx;

    // ---- toScreen: map (x, y) in data space → QPointF in widget space --
    auto toScreen = [&](double x, double y) -> QPointF {
        const double sx = pr.left + (x - xMin) / (xMax - xMin) * pr.w;
        const double sy = pr.top  + pr.h - (y - yMin) / (yMax - yMin) * pr.h;
        return {sx, sy};
    };

    // ---- Panel background -----------------------------------------------
    QRect panelRect(pr.left, pr.top, pr.w, pr.h);
    painter.fillRect(panelRect, QColor(32, 32, 38));

    // ---- Y grid lines + tick labels -------------------------------------
    const int nTicks = 4;
    {
        QFont f = painter.font();
        f.setBold(false);
        f.setPointSize(8);
        painter.setFont(f);

        for(int k = 0; k <= nTicks; ++k)
        {
            const double y = yMin + k * (yMax - yMin) / nTicks;
            const QPointF left  = toScreen(xMin, y);
            const QPointF right = toScreen(xMax, y);

            // Grid line
            painter.setPen(QPen(QColor(55, 55, 65), 1, Qt::DotLine));
            painter.drawLine(left, right);

            // Tick label
            painter.setPen(QColor(160, 160, 160));
            QString lbl;
            if(std::abs(y) < 1e-3 && std::abs(yMax - yMin) < 0.1)
                lbl = QString::number(y, 'e', 2);
            else
                lbl = QString::number(y, 'f', 3);

            painter.drawText(QRectF(0, left.y() - 8, pr.left - 4, 16),
                             Qt::AlignRight | Qt::AlignVCenter, lbl);
        }
    }

    // ---- X tick labels (on every panel for readability) ----------------
    {
        const int nXticks = 5;
        QFont f = painter.font();
        f.setPointSize(8);
        painter.setFont(f);
        painter.setPen(QColor(130, 130, 130));

        for(int k = 0; k <= nXticks; ++k)
        {
            const double x  = xMin + k * (xMax - xMin) / nXticks;
            const QPointF pt = toScreen(x, yMin);
            painter.drawText(QRectF(pt.x() - 20, pt.y() + 3, 40, 14),
                             Qt::AlignHCenter, QString::number(x, 'f', 2));
        }
    }

    // ---- Exact solution (red line) --------------------------------------
    if(!m_exactPts.empty())
    {
        painter.setPen(QPen(QColor(220, 60, 60), 1.8));
        QPointF prev;
        bool first = true;
        for(const auto& pt : m_exactPts)
        {
            QPointF sc = toScreen(pt.x, getExact(pt));
            if(!first) painter.drawLine(prev, sc);
            prev  = sc;
            first = false;
        }
    }

    // ---- FVM points (blue circles) -------------------------------------
    if(!m_fvm.empty())
    {
        painter.setPen(QPen(QColor(60, 130, 220), 1.0));
        painter.setBrush(QColor(60, 130, 220));

        for(const auto& pt : m_fvm)
        {
            QPointF sc = toScreen(pt.x, getFvm(pt));
            painter.drawEllipse(sc, 1, 1);
        }
        painter.setBrush(Qt::NoBrush);
    }

    // ---- Panel border ---------------------------------------------------
    painter.setPen(QPen(QColor(90, 90, 100), 1));
    painter.setBrush(Qt::NoBrush);
    painter.drawRect(panelRect);

    // ---- Variable label (top-left of panel) ----------------------------
    {
        QFont f = painter.font();
        f.setBold(true);
        f.setPointSize(9);
        painter.setFont(f);
        painter.setPen(QColor(200, 200, 200));
        painter.drawText(QRectF(pr.left + 6, pr.top + 4, 120, 16),
                         Qt::AlignLeft, varLabel);
    }

    // ---- Legend (first panel only) -------------------------------------
    if(varIdx == 0)
    {
        QFont f = painter.font();
        f.setBold(false);
        f.setPointSize(8);
        painter.setFont(f);

        const int lx = pr.left + pr.w - 110;
        const int ly = pr.top + 6;

        // Exact
        if(m_hasExact)
        {
            painter.setPen(QPen(QColor(220, 60, 60), 1.8));
            painter.drawLine(lx, ly + 6, lx + 20, ly + 6);
            painter.setPen(QColor(200, 200, 200));
            painter.drawText(lx + 24, ly + 10, "Exact");
        }

        // FVM
        painter.setPen(QPen(QColor(60, 130, 220), 1.0));
        painter.setBrush(QColor(60, 130, 220));
        painter.drawEllipse(QPointF(lx + 10, ly + 20), 3.0, 3.0);
        painter.setBrush(Qt::NoBrush);
        painter.setPen(QColor(200, 200, 200));
        painter.drawText(lx + 24, ly + 24, "FVM");
    }
}
