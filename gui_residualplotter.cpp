#include "gui_residualplotter.h"
#include <QPainter>
#include <cmath>
#include <algorithm>

ResidualPlotter::ResidualPlotter(QWidget* parent)
    : QWidget(parent)
{
    setMinimumSize(500, 300);
}

void ResidualPlotter::setData(const std::vector<double>& residuals)
{
    m_data = residuals;
    update();
}

void ResidualPlotter::paintEvent(QPaintEvent*)
{
    QPainter p(this);
    p.fillRect(rect(), QColor(30, 30, 30));

    if(m_data.size() < 2) return;

    const int ml = 70, mr = 20, mt = 30, mb = 40;  // margins
    QRect plot(ml, mt, width() - ml - mr, height() - mt - mb);

    // --- compute log10 min/max -----------------------------------------
    double logMin =  1e30;
    double logMax = -1e30;
    bool hasPositive = false;


    for (double v : m_data)
    {
        if (v <= 0.0)
            continue;

        hasPositive = true;
        double lv = std::log10(v);
        logMin = std::min(logMin, lv);
        logMax = std::max(logMax, lv);
    }

    if (!hasPositive)
    {
        p.setPen(QColor(180, 180, 180));
        p.drawText(rect(), Qt::AlignCenter, "No positive residuals to plot");
        return;
    }

    if (logMin >= logMax)
        logMax = logMin + 1.0;

    logMin = std::floor(logMin);
    logMax = std::ceil(logMax);

    auto toScreen = [&](int iter, double logVal) -> QPointF {
        double x = plot.left()  + (double)iter / (m_data.size()-1) * plot.width();
        double y = plot.bottom() - (logVal - logMin) / (logMax - logMin) * plot.height();
        return {x, y};
    };

    // --- grid + y labels -----------------------------------------------
    p.setPen(QColor(70, 70, 70));
    QFont f = p.font(); f.setPointSize(8); p.setFont(f);

    for(int decade = (int)logMin; decade <= (int)logMax; ++decade)
    {
        double y = plot.bottom() - (decade - logMin) / (logMax - logMin) * plot.height();
        p.setPen(QColor(70, 70, 70));
        p.drawLine(QPointF(plot.left(), y), QPointF(plot.right(), y));
        p.setPen(QColor(180, 180, 180));
        p.drawText(QRectF(0, y - 8, ml - 6, 16),
                   Qt::AlignRight | Qt::AlignVCenter,
                   QString("1e%1").arg(decade));
    }

    // --- x axis labels -------------------------------------------------
    p.setPen(QColor(180, 180, 180));
    int nTicks = 5;
    for(int t = 0; t <= nTicks; ++t)
    {
        int iter = (int)((m_data.size()-1) * t / nTicks);
        double x = plot.left() + (double)iter / (m_data.size()-1) * plot.width();
        p.drawText(QRectF(x - 25, plot.bottom() + 6, 50, 20),
                   Qt::AlignHCenter, QString::number(iter));
    }

    // --- axis titles ---------------------------------------------------
    p.save();
    p.translate(14, plot.center().y());
    p.rotate(-90);
    p.drawText(QRectF(-50, -8, 100, 16), Qt::AlignCenter, "||R||  (log10)");
    p.restore();
    p.drawText(QRectF(plot.left(), height()-18, plot.width(), 16),
               Qt::AlignCenter, "Iteration");

    // --- title ---------------------------------------------------------
    QFont fb = p.font(); fb.setBold(true); fb.setPointSize(9); p.setFont(fb);
    p.drawText(QRectF(0, 4, width(), 20), Qt::AlignCenter, "Residual History");

    // --- curve ---------------------------------------------------------
    p.setRenderHint(QPainter::Antialiasing);
    p.setPen(QPen(QColor(80, 160, 255), 1.5));

    QPointF prev;
    bool first = true;
    for(int i = 0; i < (int)m_data.size(); ++i)
    {
        if(m_data[i] <= 0) continue;
        QPointF pt = toScreen(i, std::log10(m_data[i]));
        if(!first) p.drawLine(prev, pt);
        prev  = pt;
        first = false;
    }
}
