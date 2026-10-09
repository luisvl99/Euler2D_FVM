#pragma once
#include "mesh_node.h"
#include "mesh_boundarytype.h"

struct Face
{

    int node1;
    int node2;

    int leftCell;
    int rightCell;

    double cx;
    double cy;

    double nx;
    double ny;

    double length;

    BoundaryType boundaryType;

    Face()
        :
        node1(-1),
        node2(-1),
        leftCell(-1),
        rightCell(-1),
        cx(0.0),
        cy(0.0),
        nx(0.0),
        ny(0.0),
        length(0.0),
        boundaryType(BoundaryType::Internal)
    {}

};
