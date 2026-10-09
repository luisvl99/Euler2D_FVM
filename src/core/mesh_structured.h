#pragma once

#include <vector>
#include "mesh_node.h"
#include "mesh_face.h"
#include "mesh_cell.h"

class StructuredMesh
{

public:

    int Nx;
    int Ny;

    double Lx;
    double Ly;

    std::vector<Node> nodes;
    std::vector<Cell> cells;
    std::vector<Face> faces;

    StructuredMesh(int Nx_, int Ny_, double Lx_, double Ly_);

    void generateNodes();
    void generateCells();
    void generateFaces();

    void buildCellFaceConnectivity();

    // -----------------------------------------------------------------------
    // setBoundaries
    //
    //  Retag every boundary face by side — left (nx < 0), right (nx > 0),
    //  top/bottom (nx == 0) — without rebuilding the mesh.  Call this before
    //  constructing the solver, which registers one BC per type it finds.
    // -----------------------------------------------------------------------
    void setBoundaries(BoundaryType left, BoundaryType right,
                       BoundaryType topBottom);


    int nodeIndex(int i, int j) const;
    int cellIndex(int i, int j) const;

};
