#ifndef GUI_MESHVIEWER_H
#define GUI_MESHVIEWER_H

#include <QWidget>
#include "mesh_structured.h"
#include "gui_fieldtype.h"

class MeshViewer : public QWidget
{
    Q_OBJECT

public:
    explicit MeshViewer(QWidget *parent = nullptr);

    void setMesh(StructuredMesh *mesh);

    void setDrawCellCenters(bool enabled);
    void setDrawFaceNormals(bool enabled);

    // Field visualization
    void setActiveField(FieldType ft);
    void setDrawLegend(bool enabled);

protected:
    void paintEvent(QPaintEvent *event) override;
    void wheelEvent(QWheelEvent *event) override;
    void mousePressEvent(QMouseEvent *event) override;
    void mouseMoveEvent(QMouseEvent *event) override;

private:
    void drawLegendBar(QPainter& painter, double vmin, double vmax);

    StructuredMesh *mesh = nullptr;

    bool drawCenters = false;
    bool drawNormals = false;
    bool showLegend  = true;

    FieldType activeField = FieldType::None;

    double zoom    = 1.0;
    double offsetX = 0.0;
    double offsetY = 0.0;

    QPoint lastMousePos;
};

#endif
