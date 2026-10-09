#pragma once

#include <memory>
#include <stdexcept>
#include <string>
#include <cmath>

#include "physics_eulerstate.h"
#include "physics_flowparameters.h"
#include "mesh_boundarytype.h"
#include "mesh_face.h"

// ---------------------------------------------------------------------------
// BoundaryCondition  —  abstract base class
//
//  Every concrete BC must implement exactly one method:
//
//      ghostState(U_real, face, params) → EulerState
//
//  The ghost state is a fictitious conserved state on the "outside" of the
//  boundary face.  The residual assembly loop passes (U_real, U_ghost) to
//  the numerical flux (Rusanov or HLLC) exactly as it would for an internal
//  face, so the flux's stability and conservation properties carry over.
//
//  Design rules
//  ------------
//  * Each subclass is stateless beyond what it reads from FlowParameters.
//  * No subclass modifies FlowParameters or any mesh data.
//  * Adding a new BC = adding one new subclass + one case in the factory.
//    Nothing else changes.
//
//  Normal convention
//  -----------------
//  face.nx / face.ny is the UNIT outward normal pointing FROM the interior
//  cell TOWARD the ghost cell (i.e. out of the domain).  All ghost-state
//  constructions below follow this convention.
// ---------------------------------------------------------------------------

class BoundaryCondition
{
public:

    virtual ~BoundaryCondition() = default;

    // Returns the ghost-cell conserved state for this boundary face.
    //
    // Parameters
    // ----------
    // U_real  : conserved state of the one interior cell adjacent to the face
    // face    : the boundary face (provides nx, ny and BoundaryType)
    // params  : free-stream / inlet data (read-only)
    virtual EulerState ghostState(const EulerState&    U_real,
                                  const Face&          face,
                                  const FlowParameters& params) const = 0;

    // Human-readable name — useful for debug output
    virtual const char* name() const = 0;

protected:

    // ------------------------------------------------------------------
    // Shared geometry helper: decompose velocity into normal + tangential
    //
    //  Given primitive velocity (u, v) and unit normal (nx, ny):
    //    u_n  = u·nx + v·ny          (scalar normal component)
    //    u_tx = u - u_n·nx           (tangential x-component)
    //    u_ty = v - u_n·ny           (tangential y-component)
    //
    //  This decomposition is used by Wall, Symmetry, and Farfield.
    // ------------------------------------------------------------------
    struct VelocityDecomposition
    {
        double un;   // normal component (scalar)
        double utx;  // tangential component, x
        double uty;  // tangential component, y
    };

    static VelocityDecomposition decomposeVelocity(double u, double v,
                                                   double nx, double ny)
    {
        double un  = u*nx + v*ny;
        return { un, u - un*nx, v - un*ny };
    }
};


// ===========================================================================
//  Concrete boundary condition classes
// ===========================================================================


// ---------------------------------------------------------------------------
// BC_Wall  —  inviscid (slip) wall
//
//  Physics: no flow through the wall.
//
//  Ghost state: reflect the normal velocity component, keep density,
//  tangential velocity and pressure.
//
//    u_n^ghost  = -u_n^real          (reflection)
//    u_t^ghost  =  u_t^real          (no change)
//    ρ^ghost    =  ρ^real
//    p^ghost    =  p^real
//
//  Effect on the flux: the average normal velocity at the face is
//  0.5*(u_n + (-u_n)) = 0 (Rusanov), and HLLC finds S* = 0 by symmetry, so
//  the mass flux through the wall is zero.  The pressure term still
//  contributes — this is the wall pressure force on the cell.
// ---------------------------------------------------------------------------

class BC_Wall : public BoundaryCondition
{
public:

    EulerState ghostState(const EulerState&    U_real,
                          const Face&          face,
                          const FlowParameters& /*params*/) const override
    {
        double rho = U_real.rho;
        double u   = U_real.rho_u / rho;
        double v   = U_real.rho_v / rho;
        double p   = EulerPhysics::pressure(U_real);

        auto [un, utx, uty] = decomposeVelocity(u, v, face.nx, face.ny);

        // Reflect: flip normal component, keep tangential
        double u_ghost = utx - un * face.nx;   // = u - 2*un*nx
        double v_ghost = uty - un * face.ny;   // = v - 2*un*ny

        return EulerPhysics::primitiveToConserved(rho, u_ghost, v_ghost, p);
    }

    const char* name() const override { return "Wall (slip)"; }
};


// ---------------------------------------------------------------------------
// BC_Symmetry  —  symmetry plane
//
//  Mathematically identical to the slip wall: reflect the normal velocity.
//  Kept as a separate type because it represents a different physical
//  intent (the domain is symmetric across this line, not bounded by a wall).
// ---------------------------------------------------------------------------

class BC_Symmetry : public BC_Wall
{
public:
    const char* name() const override { return "Symmetry"; }
};


// ---------------------------------------------------------------------------
// BC_Inlet  —  supersonic inlet (all characteristics enter)
//
//  Physics: at a supersonic inlet ALL four characteristics point into the
//  domain, so all four flow variables must be prescribed externally.  The
//  interior solution has no influence on the boundary state.
//
//  Ghost state: the prescribed free-stream / inlet conserved state.
//
//  Note: for a SUBSONIC inlet only three variables should be prescribed
//  (one is determined by the interior via a Riemann invariant); see
//  BC_InletSubsonic.
// ---------------------------------------------------------------------------

class BC_Inlet : public BoundaryCondition
{
public:

    EulerState ghostState(const EulerState&    /*U_real*/,
                          const Face&          /*face*/,
                          const FlowParameters& params) const override
    {
        // The interior state is deliberately ignored — the inlet fully
        // prescribes the flow.
        return params.U_inlet;
    }

    const char* name() const override { return "Inlet (supersonic)"; }
};


// ---------------------------------------------------------------------------
// BC_Outlet  —  supersonic outlet (all characteristics leave)
//
//  Physics: at a supersonic outlet ALL four characteristics point out of
//  the domain, so no information can travel back in from outside.  The
//  ghost state is simply copied from the interior cell (zero-gradient).
//
//  Ghost state: U_ghost = U_real
//
//  Note: for a SUBSONIC outlet one characteristic re-enters the domain and
//  the back-pressure must be prescribed; see BC_OutletSubsonic.
// ---------------------------------------------------------------------------

class BC_Outlet : public BoundaryCondition
{
public:

    EulerState ghostState(const EulerState&    U_real,
                          const Face&          /*face*/,
                          const FlowParameters& /*params*/) const override
    {
        return U_real;
    }

    const char* name() const override { return "Outlet (supersonic)"; }
};


// ---------------------------------------------------------------------------
// BC_Farfield  —  far-field (simple version)
//
//  Same ghost state as the supersonic inlet: the prescribed free-stream
//  state.  Exact for supersonic flow and for a uniform free stream.  A
//  characteristic (Riemann-invariant) far-field for subsonic flow is not
//  implemented; the factory would return it instead of this class.
// ---------------------------------------------------------------------------

class BC_Farfield : public BC_Inlet
{
public:
    const char* name() const override { return "Far-field (free-stream state)"; }
};



// ---------------------------------------------------------------------------
// BC_InletSubsonic  —  subsonic inlet via 1-D Riemann invariants
//
//  Characteristic analysis (outward normal n points OUT of domain,
//  flow enters so u_n = V·n < 0, |u_n| < a):
//
//    Wave speeds along n:   u_n + a > 0  (leaves domain → outgoing)
//                           u_n - a < 0  (enters domain → incoming)
//                           u_n     < 0  (enters domain → incoming)  ×2
//
//  3 incoming characteristics → prescribe from outside:
//    - Entropy:              s = p / ρ^γ  (from inlet reservoir)
//    - Tangential velocity:  u_t          (from inlet direction)
//    - Incoming invariant:   R⁻ = u_n - 2a/(γ-1)  (from inlet)
//
//  1 outgoing characteristic → comes from interior:
//    - Outgoing invariant:   R⁺ = u_n + 2a/(γ-1)  (from last cell)
//
//  Solving for the ghost state:
//    u_n^ghost = ½(R⁺_interior + R⁻_prescribed)
//    a^ghost   = ¼(γ-1)(R⁺_interior - R⁻_prescribed)
//    ρ^ghost   = (a^ghost² / (γ · s_inlet))^(1/(γ-1))
//    p^ghost   = s_inlet · ρ^ghost^γ
//    u^ghost   = u_n^ghost · nx  +  u_t (from inlet)
//    v^ghost   = u_n^ghost · ny  +  v_t (from inlet)
//
//  Note: a^ghost can become non-positive if the prescribed Mach number is
//  too high for the interior state.  A floor clamp prevents this.
// ---------------------------------------------------------------------------

class BC_InletSubsonic : public BoundaryCondition
{
public:

    EulerState ghostState(const EulerState&    U_real,
                          const Face&          face,
                          const FlowParameters& params) const override
    {
        // ------------------------------------------------------------------
        // 1. Extract interior primitive variables
        // ------------------------------------------------------------------
        const double rho_L = U_real.rho;
        const double u_L   = U_real.rho_u / rho_L;
        const double v_L   = U_real.rho_v / rho_L;
        const double a_L   = EulerPhysics::soundSpeed(U_real);
        const double un_L  = u_L * face.nx + v_L * face.ny;

        // Outgoing Riemann invariant from interior (carried by u_n + a > 0)
        const double Rplus = un_L + 2.0 * a_L / (EulerPhysics::gamma - 1.0);

        // ------------------------------------------------------------------
        // 2. Prescribed inlet state
        // ------------------------------------------------------------------
        const double rho_i = params.rho_inf;
        const double u_i   = params.u_inf;
        const double v_i   = params.v_inf;
        const double p_i   = params.p_inf;
        const double a_i   = std::sqrt(EulerPhysics::gamma * p_i / rho_i);
        const double un_i  = u_i * face.nx + v_i * face.ny;

        // Incoming Riemann invariant from the reservoir (carried by u_n - a < 0)
        const double Rminus = un_i - 2.0 * a_i / (EulerPhysics::gamma - 1.0);

        // ------------------------------------------------------------------
        // 3. Solve the 2×2 Riemann system for the ghost normal velocity
        //    and ghost sound speed
        // ------------------------------------------------------------------
        const double un_ghost = 0.5  * (Rplus + Rminus);
        double       a_ghost  = 0.25 * (EulerPhysics::gamma - 1.0) * (Rplus - Rminus);

        // Safety floor: a must be positive.  If the prescribed Mach is too
        // high, clamp to the inlet sound speed.
        if(a_ghost < 1e-6)
            a_ghost = a_i;

        // ------------------------------------------------------------------
        // 4. Ghost thermodynamic state using INLET entropy (incoming char)
        //
        //    Isentropic relation:  p = s · ρ^γ
        //    Sound speed:          a² = γ p / ρ = γ s ρ^(γ-1)
        //    →  ρ = ( a² / (γ s) )^(1/(γ-1))
        // ------------------------------------------------------------------
        const double s_inlet    = p_i / std::pow(rho_i, EulerPhysics::gamma);
        const double rho_ghost  = std::pow(
            a_ghost * a_ghost / (EulerPhysics::gamma * s_inlet),
            1.0 / (EulerPhysics::gamma - 1.0));
        const double p_ghost    = s_inlet * std::pow(rho_ghost, EulerPhysics::gamma);

        // ------------------------------------------------------------------
        // 5. Tangential velocity from the PRESCRIBED inlet (incoming char)
        // ------------------------------------------------------------------
        auto [un_unused, utx_i, uty_i] = decomposeVelocity(u_i, v_i,
                                                           face.nx, face.ny);

        // Ghost full velocity = Riemann normal + inlet tangential
        const double u_ghost = un_ghost * face.nx + utx_i;
        const double v_ghost = un_ghost * face.ny + uty_i;

        return EulerPhysics::primitiveToConserved(rho_ghost, u_ghost, v_ghost, p_ghost);
    }

    const char* name() const override { return "Inlet (subsonic, characteristic)"; }
};


// ---------------------------------------------------------------------------
// BC_OutletSubsonic  —  subsonic outlet via 1-D Riemann invariants
//
//  Characteristic analysis (outward normal n points OUT of domain,
//  flow leaves so u_n = V·n > 0, |u_n| < a):
//
//    Wave speeds along n:   u_n + a > 0  (leaves domain → outgoing)
//                           u_n     > 0  (leaves domain → outgoing)  ×2
//                           u_n - a < 0  (enters domain → incoming)
//
//  3 outgoing characteristics → come from interior:
//    - Entropy:              s = p / ρ^γ
//    - Tangential velocity:  u_t
//    - Outgoing invariant:   R⁺ = u_n + 2a/(γ-1)
//
//  1 incoming characteristic → prescribed from outside:
//    - Back pressure:  p_back  (determines R⁻ entering the domain)
//
//  Ghost state construction:
//    p^ghost   = p_back                        (prescribed)
//    ρ^ghost   = (p_back / s_interior)^(1/γ)  (interior entropy)
//    a^ghost   = sqrt(γ p_back / ρ^ghost)
//    u_n^ghost = R⁺_interior - 2 a^ghost/(γ-1)
//    u^ghost   = u_n^ghost·nx + u_t (from interior)
//    v^ghost   = u_n^ghost·ny + v_t (from interior)
//
//  Note: if p_back causes u_n^ghost to exceed a^ghost the outlet has
//  choked — physically the flow becomes sonic.  No special treatment is
//  applied here; the upstream solver will naturally settle at M ≈ 1.
// ---------------------------------------------------------------------------

class BC_OutletSubsonic : public BoundaryCondition
{
public:

    EulerState ghostState(const EulerState&    U_real,
                          const Face&          face,
                          const FlowParameters& params) const override
    {
        // ------------------------------------------------------------------
        // 1. Interior primitive variables
        // ------------------------------------------------------------------
        const double rho_L = U_real.rho;
        const double u_L   = U_real.rho_u / rho_L;
        const double v_L   = U_real.rho_v / rho_L;
        const double p_L   = EulerPhysics::pressure(U_real);
        const double a_L   = EulerPhysics::soundSpeed(U_real);
        const double un_L  = u_L * face.nx + v_L * face.ny;

        // Interior tangential velocity (outgoing characteristic)
        auto [un_unused, utx_L, uty_L] = decomposeVelocity(u_L, v_L,
                                                           face.nx, face.ny);

        // Outgoing Riemann invariant from interior
        const double Rplus = un_L + 2.0 * a_L / (EulerPhysics::gamma - 1.0);

        // ------------------------------------------------------------------
        // 2. Prescribed back-pressure → ghost thermodynamic state
        //    using INTERIOR entropy (3 outgoing chars carry entropy out)
        // ------------------------------------------------------------------
        const double p_ghost   = params.p_back;
        const double s_L       = p_L / std::pow(rho_L, EulerPhysics::gamma);
        const double rho_ghost = std::pow(p_ghost / s_L,
                                          1.0 / EulerPhysics::gamma);
        const double a_ghost   = std::sqrt(EulerPhysics::gamma * p_ghost / rho_ghost);

        // ------------------------------------------------------------------
        // 3. Ghost normal velocity from outgoing Riemann invariant
        // ------------------------------------------------------------------
        const double un_ghost = Rplus - 2.0 * a_ghost / (EulerPhysics::gamma - 1.0);

        // ------------------------------------------------------------------
        // 4. Assemble ghost velocity = Riemann normal + interior tangential
        // ------------------------------------------------------------------
        const double u_ghost = un_ghost * face.nx + utx_L;
        const double v_ghost = un_ghost * face.ny + uty_L;

        return EulerPhysics::primitiveToConserved(rho_ghost, u_ghost, v_ghost, p_ghost);
    }

    const char* name() const override { return "Outlet (subsonic, characteristic)"; }
};

// ---------------------------------------------------------------------------
// BC_ZeroGradient  —  transmissive boundary
//
//  Same ghost state as the supersonic outlet (U_ghost = U_real), used on
//  sides where waves should leave without reflection (Sod, 2-D Riemann).
// ---------------------------------------------------------------------------

class BC_ZeroGradient : public BC_Outlet
{
public:
    const char* name() const override { return "Zero-gradient (extrapolation)"; }
};

// ===========================================================================
//  BoundaryConditionFactory
//
//  Maps BoundaryType → the correct BoundaryCondition subclass.
//
//  Returns a shared_ptr so that the solver can store one instance per type
//  in a map and reuse it across all faces of that type without repeated
//  heap allocation.
//
//  To add a new BC in the future:
//    1. Write a new subclass of BoundaryCondition above.
//    2. Add one case here.
//    That is all.
// ===========================================================================

class BoundaryConditionFactory
{
public:

    static std::shared_ptr<BoundaryCondition> create(BoundaryType type)
    {
        switch(type)
        {
        case BoundaryType::Wall:
            return std::make_shared<BC_Wall>();

        case BoundaryType::Symmetry:
            return std::make_shared<BC_Symmetry>();

        case BoundaryType::Inlet:
            return std::make_shared<BC_Inlet>();

        case BoundaryType::Outlet:
            return std::make_shared<BC_Outlet>();

        case BoundaryType::Farfield:
            return std::make_shared<BC_Farfield>();

        case BoundaryType::InletSubsonic:
            return std::make_shared<BC_InletSubsonic>();

        case BoundaryType::OutletSubsonic:
            return std::make_shared<BC_OutletSubsonic>();

        case BoundaryType::ZeroGradient:
            return std::make_shared<BC_ZeroGradient>();


        case BoundaryType::Internal:
            throw std::logic_error(
                "BoundaryConditionFactory::create called with BoundaryType::Internal. "
                "Internal faces must be handled by the residual assembly loop, "
                "not by the boundary condition system.");

        default:
            throw std::logic_error(
                "BoundaryConditionFactory::create: unknown BoundaryType " +
                std::to_string(static_cast<int>(type)));
        }
    }
};
