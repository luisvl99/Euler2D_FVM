#pragma once

#include "physics_eulerstate.h"
#include "physics_primitivestate.h"

// ---------------------------------------------------------------------------
// EulerPhysics
//
//  Stateless utility class.  All methods are static.
//
//  Conventions
//  -----------
//  * The face normal (nx, ny) is assumed to be a UNIT normal pointing from
//    the left cell toward the right cell.
//  * "left" and "right" refer to the cells on each side of a face as stored
//    in Face::leftCell / Face::rightCell.
//  * The returned numerical flux F_hat has units [quantity / (m² · s)] and
//    must be multiplied by face length before accumulating into the residual.
// ---------------------------------------------------------------------------

class EulerPhysics
{
public:

    static constexpr double gamma = 1.4;

    // -----------------------------------------------------------------------
    // Basic thermodynamics
    // -----------------------------------------------------------------------

    // Convert conserved ↔ primitive variables
    static PrimitiveState conservedToPrimitive(const EulerState& U);
    static EulerState     primitiveToConserved(double rho, double u,
                                               double v,   double p);

    // Pressure from conserved state:  p = (γ-1)(ρE - ½ρ|V|²)
    static double pressure(const EulerState& U);

    // Sound speed:  a = sqrt(γ p / ρ)
    static double soundSpeed(const EulerState& U);

    // -----------------------------------------------------------------------
    // Physical flux in a given normal direction  F(U)·n
    //
    //  This is the EXACT Euler flux projected onto the face normal.
    //  It is used both directly and as a building block for numerical fluxes.
    // -----------------------------------------------------------------------
    static EulerState physicalFluxNormal(const EulerState& U,
                                         double nx,
                                         double ny);

    // -----------------------------------------------------------------------
    // Rusanov (Local Lax-Friedrichs) numerical flux
    //
    //  F_hat = 0.5 * (F(U_L,n) + F(U_R,n))  -  0.5 * S_max * (U_R - U_L)
    //
    //  where S_max = max( |u_n^L| + a^L ,  |u_n^R| + a^R )
    //
    //  Parameters
    //  ----------
    //  UL, UR : conserved states of the left and right cells
    //  nx, ny : components of the UNIT outward face normal
    //             (pointing from left cell toward right cell)
    //
    //  Returns
    //  -------
    //  The numerical interface flux F_hat (same units as physicalFluxNormal).
    //  Multiply by face length before accumulating into cell residuals.
    //
    //  Robustness guards
    //  -----------------
    //  * Density must be positive.  The solver checks every cell after each
    //    update (EulerSolver::checkStates), so a blow-up stops the run
    //    instead of reaching the fluxes.
    //  * Pressure ≤ 0 (possible in a MUSCL face state) → floored at P_MIN
    //    before the sound speed is computed.
    //  * The normal is NOT re-normalised here; the caller is responsible for
    //    supplying a unit normal (mesh_structured already does this).
    // -----------------------------------------------------------------------
    static EulerState rusanovFlux(const EulerState& UL,
                                  const EulerState& UR,
                                  double nx,
                                  double ny);

    // HLLC numerical flux (Toro §10.5).  Same conventions and guards as
    // rusanovFlux; pressures are floored at P_MIN in the wave-speed and
    // star-state formulas.
    static EulerState hllcFlux(const EulerState& UL,
                               const EulerState& UR,
                               double nx, double ny);
    // -----------------------------------------------------------------------
    // Maximum wave speed at a single state in the normal direction
    //
    //  lambda_max = |u_n| + a
    //
    //  Exposed publicly so the time-step routine can reuse it without
    //  duplicating the pressure / sound-speed logic.
    // -----------------------------------------------------------------------
    static double maxWaveSpeed(const EulerState& U,
                               double nx,
                               double ny);

private:

    // Pressure floor applied before square roots and divisions by p.
    // Intentionally small but non-zero.
    static constexpr double P_MIN = 1e-10;
};
