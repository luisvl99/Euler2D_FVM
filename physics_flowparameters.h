#pragma once

#include "physics_eulerstate.h"
#include "physics_eulerphysics.h"

// ---------------------------------------------------------------------------
// FlowParameters
//
//  Single source of truth for all externally prescribed flow data.
//  Every boundary condition class receives a const reference to one of
//  these and picks out only the fields it needs.
//
//  Fields
//  ------
//  rho_inf, u_inf, v_inf, p_inf
//      Free-stream / inlet primitive variables.
//      Used by: Inlet (supersonic), InletSubsonic (as the prescribed
//      reservoir state), Farfield.
//
//  p_back
//      Back-pressure prescribed at the subsonic outlet.
//      Used by: OutletSubsonic only.
//      For supersonic flow this field is ignored (all outlet info leaves).
//      A safe default is p_inf (no pressure jump at the exit).
//
//  U_inlet  (derived, computed once by buildDerived())
//      Conserved-variable form of the inlet state, ready for ghost-cell
//      construction in BC_Inlet and BC_Farfield without recomputing every
//      call.
//
//  Usage
//  -----
//      FlowParameters fp;
//      fp.rho_inf = 1.225;
//      fp.u_inf   = 68.0;
//      fp.v_inf   = 0.0;
//      fp.p_inf   = 101325.0;
//      fp.p_back  = 101325.0;   // equals p_inf for no pressure jump
//      fp.buildDerived();       // must be called before passing to solver
// ---------------------------------------------------------------------------

struct FlowParameters
{
    // -----------------------------------------------------------------------
    // Free-stream / inlet primitive variables (set by the user / GUI)
    // -----------------------------------------------------------------------
    double rho_inf = 1.225;      // [kg/m³]   sea-level air
    double u_inf   = 0.0;        // [m/s]     x-velocity
    double v_inf   = 0.0;        // [m/s]     y-velocity
    double p_inf   = 101325.0;   // [Pa]      sea-level pressure

    // -----------------------------------------------------------------------
    // Back-pressure for subsonic outlet  (ignored for supersonic BCs)
    // -----------------------------------------------------------------------
    double p_back  = 101325.0;   // [Pa]  default = no pressure jump at exit

    // -----------------------------------------------------------------------
    // Derived quantities (populated by buildDerived())
    // -----------------------------------------------------------------------
    EulerState U_inlet;          // conserved form of the inlet state

    // -----------------------------------------------------------------------
    // buildDerived()
    //
    //  Call this once after setting the primitive fields above.
    //  It converts primitives → conserved and stores U_inlet so that
    //  BoundaryCondition subclasses can use it directly.
    // -----------------------------------------------------------------------
    void buildDerived()
    {
        U_inlet = EulerPhysics::primitiveToConserved(rho_inf, u_inf, v_inf, p_inf);
    }
};
