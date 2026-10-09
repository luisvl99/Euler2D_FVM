#include "gui_meshviewer.h"
#include "gui_colormap.h"
#include "gui_fieldextractor.h"
#include <QPainter>
#include <QWheelEvent>
#include <QMouseEvent>
#include <QFontMetrics>

    MeshViewer::MeshViewer(QWidget *parent)
    : QWidget(parent)
{
}

void MeshViewer::setMesh(StructuredMesh *m)
{
    mesh = m;
    update();
}

void MeshViewer::setDrawCellCenters(bool enabled)
{
    drawCenters = enabled;
    update();
}

void MeshViewer::setDrawFaceNormals(bool enabled)
{
    drawNormals = enabled;
    update();
}

void MeshViewer::setActiveField(FieldType ft)
{
    activeField = ft;
    update();
}

void MeshViewer::setDrawLegend(bool enabled)
{
    showLegend = enabled;
    update();
}

void MeshViewer::paintEvent(QPaintEvent *)
{
    QPainter painter(this);
    painter.fillRect(rect(), Qt::white);

    if(!mesh)
        return;

    double Lx = mesh->Lx;
    double Ly = mesh->Ly;

    double scaleX = width()  / Lx;
    double scaleY = height() / Ly;
    double scaleFactor = std::min(scaleX, scaleY);

    double cx = Lx * 0.5;
    double cy = Ly * 0.5;
    double screenCx = width()  * 0.5;
    double screenCy = height() * 0.5;

    // -----------------------------------------------------------------------
    // Pre-compute field range (before any painter transform so we can also
    // pass it to the legend which is drawn in screen coordinates)
    // -----------------------------------------------------------------------
    FieldExtractor::Range fieldRange = {0.0, 1.0};
    bool coloringEnabled = (activeField != FieldType::None);

    if(coloringEnabled)
        fieldRange = FieldExtractor::computeRange(mesh->cells, activeField);

    // -----------------------------------------------------------------------
    // World transform
    // -----------------------------------------------------------------------
    painter.save();
    painter.translate(screenCx + offsetX, screenCy + offsetY);
    painter.scale(scaleFactor * zoom, -scaleFactor * zoom);
    painter.translate(-cx, -cy);

    // -----------------------------------------------------------------------
    // Draw cells (filled if a field is active, otherwise white)
    // -----------------------------------------------------------------------
    painter.setRenderHint(QPainter::Antialiasing, false);

    for(const Cell& cell : mesh->cells)
    {
        QPolygonF poly;
        for(int nodeID : cell.nodes)
        {
            const Node& node = mesh->nodes[nodeID];
            poly << QPointF(node.x, node.y);
        }

        if(coloringEnabled)
        {
            double val = FieldExtractor::extract(cell, activeField);
            double t   = ColorMap::normalize(val, fieldRange.vmin, fieldRange.vmax);
            QColor c   = FieldExtractor::isDiverging(activeField)
                           ? ColorMap::diverging(t)
                           : ColorMap::sequential(t);

            painter.setPen(Qt::NoPen);        // hide edges when colored
            painter.setBrush(c);
        }
        else
        {
            painter.setPen(QPen(Qt::black, 0));
            painter.setBrush(Qt::NoBrush);
        }

        painter.drawPolygon(poly);
    }

    // Thin black grid on top when field is active (helps keep orientation)
    if(coloringEnabled)
    {
        painter.setBrush(Qt::NoBrush);
        painter.setPen(QPen(QColor(0,0,0,40), 0));   // very subtle

        for(const Cell& cell : mesh->cells)
        {
            QPolygonF poly;
            for(int nodeID : cell.nodes)
            {
                const Node& nd = mesh->nodes[nodeID];
                poly << QPointF(nd.x, nd.y);
            }
            painter.drawPolygon(poly);
        }
    }

    // -----------------------------------------------------------------------
    // Cell centers
    // -----------------------------------------------------------------------
    if(drawCenters)
    {
        painter.setBrush(Qt::red);
        painter.setPen(Qt::NoPen);
        for(const Cell& cell : mesh->cells)
            painter.drawEllipse(QPointF(cell.cx, cell.cy), 0.01, 0.01);
    }

    // -----------------------------------------------------------------------
    // Face normals
    // -----------------------------------------------------------------------
    if(drawNormals)
    {
        double cellSize = std::min(mesh->Lx / mesh->Nx,
                                   mesh->Ly / mesh->Ny);
        double scale = cellSize * 0.35;
        painter.setPen(QPen(Qt::blue, 0));

        for(const Face& face : mesh->faces)
        {
            double startX = face.cx - 0.5 * scale * face.nx;
            double startY = face.cy - 0.5 * scale * face.ny;
            double endX   = face.cx + 0.5 * scale * face.nx;
            double endY   = face.cy + 0.5 * scale * face.ny;

            painter.drawLine(QPointF(startX, startY), QPointF(endX, endY));
            painter.setBrush(Qt::blue);
            painter.drawEllipse(QPointF(endX, endY), scale*0.08, scale*0.08);
            painter.setBrush(Qt::NoBrush);
        }
    }

    painter.restore();

    // -----------------------------------------------------------------------
    // Legend (screen coordinates — drawn AFTER restoring transform)
    // -----------------------------------------------------------------------
    if(coloringEnabled && showLegend)
        drawLegendBar(painter, fieldRange.vmin, fieldRange.vmax);
}

// ---------------------------------------------------------------------------
// Draw a vertical color-bar legend in the bottom-right corner
// ---------------------------------------------------------------------------
void MeshViewer::drawLegendBar(QPainter& painter, double vmin, double vmax)
{
    const int barW  = 18;
    const int barH  = 160;
    const int margR = 16;   // right margin
    const int margB = 16;   // bottom margin
    const int lblW  = 60;   // space reserved for text left of bar

    int x0 = width()  - margR - barW - lblW;
    int y0 = height() - margB - barH;

    // Background panel
    QRect panel(x0 - 10, y0 - 24, barW + lblW + 20, barH + 38);
    painter.setPen(QPen(QColor(180,180,180), 1));
    painter.setBrush(QColor(255,255,255,210));
    painter.drawRoundedRect(panel, 6, 6);

    // Color bar — draw as thin horizontal strips
    int N = barH;
    for(int i = 0; i < N; i++)
    {
        double t = 1.0 - (double)i / (N - 1);   // top = max, bottom = min
        QColor c = FieldExtractor::isDiverging(activeField)
                       ? ColorMap::diverging(t)
                       : ColorMap::sequential(t);
        painter.setPen(QPen(c, 1));
        painter.drawLine(x0 + lblW,     y0 + i,
                         x0 + lblW + barW - 1, y0 + i);
    }

    // Border around bar
    painter.setPen(QPen(Qt::gray, 1));
    painter.setBrush(Qt::NoBrush);
    painter.drawRect(x0 + lblW, y0, barW, barH);

    // Tick labels: max (top), mid, min (bottom)
    painter.setPen(Qt::black);
    QFont f = painter.font();
    f.setPointSize(8);
    painter.setFont(f);
    QFontMetrics fm(f);

    auto fmtVal = [](double v) -> QString {
        if(std::abs(v) < 1e-3 || std::abs(v) >= 1e5)
            return QString::number(v, 'e', 2);
        return QString::number(v, 'f', 4);
    };

    QString smax = fmtVal(vmax);
    QString smid = fmtVal(0.5*(vmin+vmax));
    QString smin = fmtVal(vmin);

    int tx = x0 + lblW - 4;
    painter.drawText(tx - fm.horizontalAdvance(smax), y0 + fm.ascent(),                    smax);
    painter.drawText(tx - fm.horizontalAdvance(smid), y0 + barH/2 + fm.ascent()/2,         smid);
    painter.drawText(tx - fm.horizontalAdvance(smin), y0 + barH - 2,                        smin);

    // Field name label above bar
    static const QMap<int,QString> fieldNames = {
                                                  {(int)FieldType::Density,     "ρ (density)"},
                                                  {(int)FieldType::VelocityMag, "|V| (m/s)"},
                                                  {(int)FieldType::VelocityU,   "u (m/s)"},
                                                  {(int)FieldType::VelocityV,   "v (m/s)"},
                                                  {(int)FieldType::Pressure,    "p (Pa)"},
                                                  {(int)FieldType::MachNumber,  "Mach"},
                                                  {(int)FieldType::TotalEnergy, "E (J/kg)"},
                                                  {(int)FieldType::Residual,    "|R|"},
                                                  };

    QString name = fieldNames.value((int)activeField, "field");
    QFont fb = painter.font();
    fb.setBold(true);
    painter.setFont(fb);
    int nx = x0 + lblW + barW/2 - fm.horizontalAdvance(name)/2;
    painter.drawText(nx, y0 - 8, name);
}

// ---------------------------------------------------------------------------
// Mouse / wheel
// ---------------------------------------------------------------------------
void MeshViewer::wheelEvent(QWheelEvent *event)
{
    double factor = 1.15;
    if(event->angleDelta().y() > 0) zoom *= factor;
    else                             zoom /= factor;
    update();
}

void MeshViewer::mousePressEvent(QMouseEvent *event)
{
    lastMousePos = event->pos();
}

void MeshViewer::mouseMoveEvent(QMouseEvent *event)
{
    if (!(event->buttons() & Qt::LeftButton)) return;
    QPoint delta = event->pos() - lastMousePos;
    offsetX += delta.x();
    offsetY += delta.y();
    lastMousePos = event->pos();
    update();
}
