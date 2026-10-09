#pragma once
#include <vector>
#include "physics_eulerstate.h"

struct Cell
{

    std::vector<int> nodes;
    std::vector<int> faces;

    int i;
    int j;

    double cx;
    double cy;

    double area;

    EulerState U;
    EulerState R;

    Cell()
        :
        i(0),
        j(0),
        cx(0.0),
        cy(0.0),
        area(0.0)
    {}

};
