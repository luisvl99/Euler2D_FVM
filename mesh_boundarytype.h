#pragma once
enum class BoundaryType
{
    Internal,
    Wall,
    Inlet,
    Outlet,
    Symmetry,
    Farfield,
    InletSubsonic,
    OutletSubsonic,
    ZeroGradient
};
