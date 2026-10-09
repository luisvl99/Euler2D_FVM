#pragma once

enum class FieldType
{
    None,           // wireframe only (default)
    Density,        // rho
    VelocityMag,    // sqrt(u^2 + v^2)
    VelocityU,      // u  (can be negative → diverging colormap)
    VelocityV,      // v
    Pressure,       // p
    MachNumber,     // |V| / a
    TotalEnergy,    // rho_E / rho  (specific total energy)
    Residual        // |R| magnitude — useful once the solver runs
};
