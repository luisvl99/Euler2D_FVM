#include "physics_eulerphysics.h"
#include <algorithm>   // std::max
#include <cmath>       // std::sqrt, std::abs

// ===========================================================================
//  Basic thermodynamics
// ===========================================================================

PrimitiveState EulerPhysics::conservedToPrimitive(const EulerState& U)
{
    PrimitiveState W;

    W.rho = U.rho;
    W.u   = U.rho_u / U.rho;
    W.v   = U.rho_v / U.rho;

    double kinetic = 0.5 * U.rho * (W.u*W.u + W.v*W.v);

    W.p = (gamma - 1.0) * (U.rho_E - kinetic);

    return W;
}

EulerState EulerPhysics::primitiveToConserved(double rho, double u,
                                              double v,   double p)
{
    return { rho, rho * u, rho * v,
             p / (gamma - 1.0) + 0.5 * rho * (u*u + v*v) };
}

double EulerPhysics::pressure(const EulerState& U)
{
    double u = U.rho_u / U.rho;
    double v = U.rho_v / U.rho;

    double kinetic = 0.5 * U.rho * (u*u + v*v);

    return (gamma - 1.0) * (U.rho_E - kinetic);
}

double EulerPhysics::soundSpeed(const EulerState& U)
{
    double p = pressure(U);

    // Guard: clamp pressure before sqrt to avoid NaN
    p = std::max(p, P_MIN);

    return std::sqrt(gamma * p / U.rho);
}

// ===========================================================================
//  Physical flux projected onto face normal:  F(U) · n
//
//  For the 2-D Euler equations the flux vector in the normal direction is:
//
//    F · n = [ ρ u_n                          ]
//            [ ρ u u_n  +  p nx               ]
//            [ ρ v u_n  +  p ny               ]
//            [ (ρE + p) u_n                   ]
//
//  where u_n = u nx + v ny  is the normal velocity component.
// ===========================================================================

EulerState EulerPhysics::physicalFluxNormal(const EulerState& U,
                                            double nx,
                                            double ny)
{
    PrimitiveState W = conservedToPrimitive(U);

    double un = W.u*nx + W.v*ny;
    double p  = W.p;
    EulerState F;
    F.rho   = U.rho   * un;
    F.rho_u = U.rho_u * un + p * nx;
    F.rho_v = U.rho_v * un + p * ny;
    F.rho_E = (U.rho_E + p) * un;

    return F;
}

// ===========================================================================
//  Maximum wave speed at a single state in direction (nx, ny)
//
//    λ_max = |u_n| + a
//
//  This is the largest eigenvalue magnitude of the Euler Jacobian projected
//  onto the face normal.  It appears in both the Rusanov flux and the CFL
//  time-step restriction.
// ===========================================================================

double EulerPhysics::maxWaveSpeed(const EulerState& U,
                                  double nx,
                                  double ny)
{
    double u  = U.rho_u / U.rho;
    double v  = U.rho_v / U.rho;
    double un = u*nx + v*ny;
    double a  = soundSpeed(U);   // already guards negative pressure

    return std::abs(un) + a;
}

// ===========================================================================
//  Rusanov (Local Lax-Friedrichs) numerical flux
//
//  Formula:
//
//    F_hat = 0.5 * ( F(U_L, n) + F(U_R, n) )
//          - 0.5 * S_max * ( U_R - U_L )
//
//  where
//
//    S_max = max( |u_n^L| + a^L ,  |u_n^R| + a^R )
//
//  Derivation sketch
//  -----------------
//  Consider a 1-D Riemann problem locally at the face.  The Rusanov scheme
//  bounds the fastest signal speed from EITHER side by S_max and adds a
//  scalar dissipation proportional to that speed times the state jump.
//  This guarantees:
//    (a) Conservation  – flux is single-valued at the interface
//    (b) Consistency   – reduces to F(U,n) when U_L = U_R
//    (c) Entropy stability – the added diffusion prevents expansion shocks
//
//  Sign convention
//  ---------------
//  The returned flux is defined as positive in the direction of (nx, ny).
//  When accumulating into cell residuals:
//    * LEFT  cell:  residual -= F_hat * face.length   (flux leaves left cell)
//    * RIGHT cell:  residual += F_hat * face.length   (flux enters right cell)
//  (The residual assembly loop in the solver handles this.)
// ===========================================================================

EulerState EulerPhysics::rusanovFlux(const EulerState& UL,
                                     const EulerState& UR,
                                     double nx,
                                     double ny)
{
    // ------------------------------------------------------------------
    // 1. Physical fluxes on each side:  F(U_L)·n  and  F(U_R)·n
    // ------------------------------------------------------------------
    EulerState FL = physicalFluxNormal(UL, nx, ny);
    EulerState FR = physicalFluxNormal(UR, nx, ny);

    // ------------------------------------------------------------------
    // 2. Maximum wave speed (local, evaluated independently on each side)
    //
    //    S_L = |u_n^L| + a^L
    //    S_R = |u_n^R| + a^R
    //    S_max = max(S_L, S_R)
    // ------------------------------------------------------------------
    double SL   = maxWaveSpeed(UL, nx, ny);
    double SR   = maxWaveSpeed(UR, nx, ny);
    double Smax = std::max(SL, SR);

    // ------------------------------------------------------------------
    // 3. Rusanov flux
    //
    //    F_hat = 0.5*(FL + FR)  -  0.5*Smax*(UR - UL)
    //          = central part   -  dissipation part
    // ------------------------------------------------------------------
    EulerState Fhat = (FL + FR) * 0.5  -  (UR - UL) * (0.5 * Smax);

    return Fhat;
}

EulerState EulerPhysics::hllcFlux(const EulerState& UL,
                                  const EulerState& UR,
                                  double nx, double ny)
{
    // ---- Primitives -------------------------------------------------------
    // Pressures are floored at P_MIN, as in soundSpeed(): a MUSCL face state
    // can have p <= 0, and pK is divided by in qFactor and in the star state.
    const double rhoL = UL.rho,  rhoR = UR.rho;
    const double uL   = UL.rho_u / rhoL, uR = UR.rho_u / rhoR;
    const double vL   = UL.rho_v / rhoL, vR = UR.rho_v / rhoR;
    const double pL   = std::max(pressure(UL), P_MIN);
    const double pR   = std::max(pressure(UR), P_MIN);
    const double aL   = soundSpeed(UL),  aR = soundSpeed(UR);
    const double EL   = UL.rho_E / rhoL, ER = UR.rho_E / rhoR; // specific total energy

    // Normal velocities
    const double unL = uL*nx + vL*ny;
    const double unR = uR*nx + vR*ny;

    // ---- Wave speed estimates (pressure-based, Toro §10.5) ----------------
    // PVRS pressure estimate used to compute left/right factors qL, qR
    const double rhoBar = 0.5*(rhoL + rhoR);
    const double aBar   = 0.5*(aL   + aR);
    const double pPVRS  = std::max(0.0,
                                  0.5*(pL + pR) - 0.5*(unR - unL)*rhoBar*aBar);

    auto qFactor = [&](double p_star, double pK) {
        if(p_star <= pK) return 1.0;
        return std::sqrt(1.0 + (gamma+1.0)/(2.0*gamma) * (p_star/pK - 1.0));
    };

    const double SL = unL - aL * qFactor(pPVRS, pL);
    const double SR = unR + aR * qFactor(pPVRS, pR);

    // ---- Contact wave speed S* --------------------------------------------
    const double num   = pR - pL + rhoL*unL*(SL - unL) - rhoR*unR*(SR - unR);
    const double denom = rhoL*(SL - unL) - rhoR*(SR - unR);
    const double Sstar = (std::abs(denom) > 1e-30) ? num/denom : 0.5*(unL+unR);

    // ---- Physical fluxes --------------------------------------------------
    const EulerState FL = physicalFluxNormal(UL, nx, ny);
    const EulerState FR = physicalFluxNormal(UR, nx, ny);

    // ---- Helper: build U*K star state -------------------------------------
    auto starState = [&](double rhoK, double unK,
                         double uK,   double vK,
                         double pK,   double EK,
                         double SK) -> EulerState
    {
        const double factor = rhoK * (SK - unK) / (SK - Sstar);
        EulerState Ustar;
        Ustar.rho   = factor;
        Ustar.rho_u = factor * (uK + (Sstar - unK)*nx);
        Ustar.rho_v = factor * (vK + (Sstar - unK)*ny);
        Ustar.rho_E = factor * (EK + (Sstar - unK)*
                                         (Sstar + pK / (rhoK*(SK - unK))));
        return Ustar;
    };

    // ---- Select region and return flux ------------------------------------
    if(SL >= 0.0)
    {
        return FL;
    }
    else if(Sstar >= 0.0)   // SL < 0 <= S*  → left star region
    {
        EulerState UstarL = starState(rhoL, unL, uL, vL, pL, EL, SL);
        return FL + (UstarL - UL) * SL;
    }
    else if(SR >= 0.0)      // S* < 0 <= SR  → right star region
    {
        EulerState UstarR = starState(rhoR, unR, uR, vR, pR, ER, SR);
        return FR + (UstarR - UR) * SR;
    }
    else
    {
        return FR;
    }
}
