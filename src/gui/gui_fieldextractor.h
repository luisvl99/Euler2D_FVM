#pragma once
#include "gui_fieldtype.h"
#include "mesh_cell.h"
#include "physics_eulerphysics.h"
#include <cmath>
#include <limits>
#include <vector>

// ---------------------------------------------------------------------------
// FieldExtractor
//
//  Given a FieldType, extracts a scalar value from a single cell and also
//  computes the global [min, max] range across the whole mesh in one pass.
//
//  FieldType::VelocityU / VelocityV use diverging colormap (signed);
//  all others use sequential.
// ---------------------------------------------------------------------------

class FieldExtractor
{
public:

    static double extract(const Cell& cell, FieldType ft)
    {
        const EulerState& U = cell.U;

        if(U.rho < 1e-30)
            return 0.0;

        switch(ft)
        {
        case FieldType::Density:
            return U.rho;

        case FieldType::VelocityMag:
        {
            double u = U.rho_u / U.rho;
            double v = U.rho_v / U.rho;
            return std::sqrt(u*u + v*v);
        }

        case FieldType::VelocityU:
            return U.rho_u / U.rho;

        case FieldType::VelocityV:
            return U.rho_v / U.rho;

        case FieldType::Pressure:
            return EulerPhysics::pressure(U);

        case FieldType::MachNumber:
        {
            double u = U.rho_u / U.rho;
            double v = U.rho_v / U.rho;
            double spd = std::sqrt(u*u + v*v);
            double a   = EulerPhysics::soundSpeed(U);
            return (a > 1e-30) ? spd / a : 0.0;
        }

        case FieldType::TotalEnergy:
            return U.rho_E / U.rho;

        case FieldType::Residual:
        {
            // |R| = sqrt(R0^2 + R1^2 + R2^2 + R3^2)
            const EulerState& R = cell.R;
            return std::sqrt(R.rho*R.rho
                             + R.rho_u*R.rho_u
                             + R.rho_v*R.rho_v
                             + R.rho_E*R.rho_E);
        }

        default:
            return 0.0;
        }
    }

    // Compute global min/max for the field across all cells
    struct Range { double vmin; double vmax; };

    static Range computeRange(const std::vector<Cell>& cells, FieldType ft)
    {
        double vmin =  std::numeric_limits<double>::max();
        double vmax = -std::numeric_limits<double>::max();

        for(const Cell& c : cells)
        {
            double v = extract(c, ft);
            if(v < vmin) vmin = v;
            if(v > vmax) vmax = v;
        }

        // Avoid degenerate range (e.g. uniform field)
        if(std::abs(vmax - vmin) < 1e-30)
        {
            vmin -= 0.5;
            vmax += 0.5;
        }

        return {vmin, vmax};
    }

    // True for fields that should use the diverging (blue-white-red) colormap
    static bool isDiverging(FieldType ft)
    {
        return ft == FieldType::VelocityU
               || ft == FieldType::VelocityV;
    }
};
