#include "mesh_structured.h"
#include "mesh_boundarytype.h"
#include <cmath>

StructuredMesh::StructuredMesh(int Nx_, int Ny_, double Lx_, double Ly_)
{
    Nx = Nx_;
    Ny = Ny_;
    Lx = Lx_;
    Ly = Ly_;

    generateNodes();
    generateCells();
    generateFaces();

    buildCellFaceConnectivity();
}

void StructuredMesh::generateNodes()
{

    double dx = Lx / Nx;
    double dy = Ly / Ny;

    for(int j=0; j<=Ny; j++)
    {
        for(int i=0; i<=Nx; i++)
        {
            nodes.emplace_back(i*dx, j*dy);
        }
    }

}

void StructuredMesh::generateCells()
{
    for(int j=0; j<Ny; j++)
    {
        for(int i=0; i<Nx; i++)
        {

            Cell cell;

            cell.i = i;
            cell.j = j;

            int n1 = nodeIndex(i,j);
            int n2 = nodeIndex(i+1,j);
            int n3 = nodeIndex(i+1,j+1);
            int n4 = nodeIndex(i,j+1);

            cell.nodes = {n1,n2,n3,n4};

            Node &A = nodes[n1];
            Node &C = nodes[n3];

            cell.cx = 0.5*(A.x + C.x);
            cell.cy = 0.5*(A.y + C.y);

            double dx = nodes[n2].x - nodes[n1].x;
            double dy = nodes[n4].y - nodes[n1].y;

            cell.area = dx*dy;

            cells.push_back(cell);


        }
    }
}

void StructuredMesh::generateFaces()
{


    // vertical faces
    for(int j=0; j<Ny; j++)
    {
        for(int i=0; i<=Nx; i++)
        {

            Face f;

            int n1 = nodeIndex(i,j);
            int n2 = nodeIndex(i,j+1);

            f.node1 = n1;
            f.node2 = n2;

            Node &A = nodes[n1];
            Node &B = nodes[n2];

            f.cx = 0.5*(A.x + B.x);
            f.cy = 0.5*(A.y + B.y);

            double dx = B.x - A.x;
            double dy = B.y - A.y;

            f.length = std::sqrt(dx*dx + dy*dy);

            f.nx = 1.0;
            f.ny = 0.0;

            // Default tags (inlet left, outlet right, walls top/bottom);
            // test cases retag them with setBoundaries()
            if(i == 0)
            {
                // Left boundary: outward normal points LEFT
                f.nx = -1.0;
                f.boundaryType = BoundaryType::Inlet;
            }
            else if(i == Nx)
            {
                // Right boundary: outward normal points RIGHT
                f.boundaryType = BoundaryType::Outlet;
            }

            if(i>0)
                f.leftCell = cellIndex(i-1,j);

            if(i<Nx)
                f.rightCell = cellIndex(i,j);

            faces.push_back(f);

        }
    }

    // horizontal faces
    for(int j=0; j<=Ny; j++)
    {
        for(int i=0; i<Nx; i++)
        {
            Face f;

            int n1 = nodeIndex(i,   j);
            int n2 = nodeIndex(i+1, j);

            f.node1 = n1;
            f.node2 = n2;

            Node &A = nodes[n1];
            Node &B = nodes[n2];

            f.cx = 0.5*(A.x + B.x);
            f.cy = 0.5*(A.y + B.y);

            double ddx = B.x - A.x;
            double ddy = B.y - A.y;

            f.length = std::sqrt(ddx*ddx + ddy*ddy);

            f.nx = 0.0;
            f.ny = 1.0;

            if(j == 0)
            {
                // Bottom boundary: outward normal points DOWN
                f.ny = -1.0;
                f.boundaryType = BoundaryType::Wall;
            }
            else if(j == Ny)
            {
                // Top boundary: outward normal points UP
                f.boundaryType = BoundaryType::Wall;
            }

            if(j > 0)
                f.leftCell  = cellIndex(i, j-1);

            if(j < Ny)
                f.rightCell = cellIndex(i, j);

            faces.push_back(f);
        }
    }

}

void StructuredMesh::buildCellFaceConnectivity()
{

    for(int f = 0; f < static_cast<int>(faces.size()); ++f)
    {

        int L = faces[f].leftCell;
        int R = faces[f].rightCell;

        if(L >= 0)
            cells[L].faces.push_back(f);

        if(R >= 0)
            cells[R].faces.push_back(f);

    }

}


// ---------------------------------------------------------------------------
// setBoundaries
//
//  The side is taken from the face normal rather than the current tag, so
//  any sequence of test cases can be run on the same mesh.
// ---------------------------------------------------------------------------
void StructuredMesh::setBoundaries(BoundaryType left, BoundaryType right,
                                   BoundaryType topBottom)
{
    for(Face& f : faces)
    {
        if(f.boundaryType == BoundaryType::Internal) continue;

        if(f.nx < 0.0)      f.boundaryType = left;
        else if(f.nx > 0.0) f.boundaryType = right;
        else                f.boundaryType = topBottom;
    }
}

int StructuredMesh::nodeIndex(int i, int j) const
{
    return j*(Nx+1) + i;
}

int StructuredMesh::cellIndex(int i, int j) const
{
    return j*Nx + i;
}



