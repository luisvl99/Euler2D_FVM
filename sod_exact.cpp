#include "sod_exact.h"
#include <algorithm>
#include <iostream>

// ===========================================================================
//  Constructor / setStates
// ===========================================================================

SodExact::SodExact(State L, State R, double gamma)
{
    setStates(L, R, gamma);
}

void SodExact::setStates(State L, State R, double gamma)
{
    m_L      = L;
    m_R      = R;
    m_gamma  = gamma;
    m_solved = false;
    solve();
}

// ===========================================================================
//  Pressure function  f(p, K)
//
//  This encodes how much velocity change a wave of pressure p produces on
//  state K.  F(p) = f(p,L) + f(p,R) + (uR - uL) = 0 is the equation whose
//  root gives p*.
//
//  Shock   (p > pK):
//    f = (p - pK) * sqrt(AK / (p + BK))
//    AK = 2 / ((γ+1) ρK),   BK = (γ-1)/(γ+1) * pK
//
//  Rarefaction  (p ≤ pK):
//    f = 2aK/(γ-1) * ((p/pK)^((γ-1)/(2γ)) - 1)
// ===========================================================================

double SodExact::f(double p, const State& K, double aK) const
{
    const double g = m_gamma;
    if(p > K.p)
    {
        // Shock relation
        const double A = 2.0 / ((g + 1.0) * K.rho);
        const double B = (g - 1.0) / (g + 1.0) * K.p;
        return (p - K.p) * std::sqrt(A / (p + B));
    }
    else
    {
        // Isentropic rarefaction
        return (2.0 * aK / (g - 1.0))
               * (std::pow(p / K.p, (g - 1.0) / (2.0 * g)) - 1.0);
    }
}

// ===========================================================================
//  df/dp  —  derivative of f for Newton's method
//
//  Shock   (p > pK):
//    df/dp = sqrt(AK/(p+BK)) * (1 - (p-pK) / (2*(p+BK)))
//
//  Rarefaction  (p ≤ pK):
//    df/dp = 1/(ρK aK) * (p/pK)^(-(γ+1)/(2γ))
// ===========================================================================

double SodExact::dfdp(double p, const State& K, double aK) const
{
    const double g = m_gamma;
    if(p > K.p)
    {
        const double A = 2.0 / ((g + 1.0) * K.rho);
        const double B = (g - 1.0) / (g + 1.0) * K.p;
        return std::sqrt(A / (p + B)) * (1.0 - (p - K.p) / (2.0 * (p + B)));
    }
    else
    {
        return (1.0 / (K.rho * aK))
        * std::pow(p / K.p, -(g + 1.0) / (2.0 * g));
    }
}

// ===========================================================================
//  solve()  —  Newton iteration to find p* and u*
//
//  Initial guess: Toro's PVRS approximation (Eq. 9.28 in Toro 3rd ed.)
//  which is a linearisation around the average sound speed.  Works well
//  for weak to moderate jumps; Newton converges rapidly from there.
// ===========================================================================

void SodExact::solve()
{
    const double g  = m_gamma;
    const double aL = soundSpeed(m_L, g);
    const double aR = soundSpeed(m_R, g);

    // PVRS initial guess
    double p = 0.5 * (m_L.p + m_R.p)
               - 0.125 * (m_R.u - m_L.u) * (m_L.rho + m_R.rho) * (aL + aR);
    p = std::max(p, 1e-10);   // pressure floor

    const int    MAXITER = 100;
    const double TOL     = 1e-10;

    bool converged = false;
    for(int iter = 0; iter < MAXITER; ++iter)
    {
        const double FL  = f(p, m_L, aL);
        const double FR  = f(p, m_R, aR);
        const double F   = FL + FR + (m_R.u - m_L.u);
        const double dF  = dfdp(p, m_L, aL) + dfdp(p, m_R, aR);

        if(std::abs(dF) < 1e-30) break;

        const double dp  = -F / dF;
        const double p_new = std::max(p + dp, 1e-10);

        // Relative convergence test
        if(std::abs(p_new - p) / (0.5 * (p_new + p)) < TOL)
        {
            p = p_new;
            converged = true;
            break;
        }
        p = p_new;
    }

    if(!converged)
        std::cerr << "[SodExact] WARNING: Newton iteration did not converge."
                  << "  p* = " << p << "\n";

    m_pStar = p;
    m_uStar = 0.5 * (m_L.u + m_R.u)
              + 0.5 * (f(p, m_R, aR) - f(p, m_L, aL));
    m_solved = true;

    std::cout << "[SodExact] p* = " << m_pStar
              << "  u* = " << m_uStar << "\n";
}

// ===========================================================================
//  starDensity
//
//  Density in the star region on the side of state K.
//
//  Shock  (p* > pK)  — Rankine-Hugoniot:
//    ρ* = ρK * (p*/pK + (γ-1)/(γ+1)) / ((γ-1)/(γ+1) * p*/pK + 1)
//
//  Rarefaction  (p* ≤ pK)  — isentropic:
//    ρ* = ρK * (p*/pK)^(1/γ)
// ===========================================================================

double SodExact::starDensity(double pStar, const State& K) const
{
    const double g = m_gamma;
    if(pStar > K.p)
    {
        const double ratio = pStar / K.p;
        const double mu    = (g - 1.0) / (g + 1.0);
        return K.rho * (ratio + mu) / (mu * ratio + 1.0);
    }
    else
    {
        return K.rho * std::pow(pStar / K.p, 1.0 / g);
    }
}

// ===========================================================================
//  sample()
//
//  Evaluate the exact solution at position x, time t, with the initial
//  discontinuity located at x0.
//
//  Self-similar variable:  ξ = (x − x0) / t
//
//  Wave structure (left → right):
//
//    Region 1 | Left wave (shock or fan) | 3* | contact | 4* | Right wave | 5
//
//  Left-wave bounds [SL_left, SL_right]:
//    shock      → both equal the shock speed  SL
//    rarefaction → SL_left  = uL − aL       (head, leftmost)
//                  SL_right = u* − a*L      (tail, rightmost)
//
//  Right-wave bounds [SR_left, SR_right]:
//    shock      → both equal the shock speed  SR
//    rarefaction → SR_left  = u* + a*R      (tail, leftmost)
//                  SR_right = uR + aR       (head, rightmost)
// ===========================================================================

SodExact::Point SodExact::sample(double x, double t, double x0) const
{
    // Edge case: t = 0 → return initial condition
    if(t < 1e-14)
    {
        if(x < x0) return {x, m_L.rho, m_L.u, m_L.p};
        else        return {x, m_R.rho, m_R.u, m_R.p};
    }

    const double g     = m_gamma;
    const double aL    = soundSpeed(m_L, g);
    const double aR    = soundSpeed(m_R, g);
    const double xi    = (x - x0) / t;

    // Star-region densities and sound speeds
    const double rhoSL = starDensity(m_pStar, m_L);
    const double rhoSR = starDensity(m_pStar, m_R);
    const double aSL   = std::sqrt(g * m_pStar / rhoSL);
    const double aSR   = std::sqrt(g * m_pStar / rhoSR);

    // ---- Left-wave bounds ----
    double SL_left, SL_right;
    if(m_pStar > m_L.p)
    {
        // Left shock speed (Eq. 10.54 in Toro)
        const double SL = m_L.u - aL * std::sqrt((g + 1.0) / (2.0 * g)
                                                     * m_pStar / m_L.p
                                                 + (g - 1.0) / (2.0 * g));
        SL_left = SL_right = SL;
    }
    else
    {
        // Left rarefaction
        SL_left  = m_L.u - aL;      // head (fastest leftward)
        SL_right = m_uStar - aSL;   // tail (slowest leftward)
    }

    // ---- Right-wave bounds ----
    double SR_left, SR_right;
    if(m_pStar > m_R.p)
    {
        // Right shock speed
        const double SR = m_R.u + aR * std::sqrt((g + 1.0) / (2.0 * g)
                                                     * m_pStar / m_R.p
                                                 + (g - 1.0) / (2.0 * g));
        SR_left = SR_right = SR;
    }
    else
    {
        // Right rarefaction
        SR_left  = m_uStar + aSR;   // tail (slowest rightward)
        SR_right = m_R.u + aR;      // head (fastest rightward)
    }

    // ---- Determine region and return primitive state ----

    if(xi <= SL_left)
    {
        // Region 1 — undisturbed left state
        return {x, m_L.rho, m_L.u, m_L.p};
    }
    else if(xi <= SL_right)
    {
        // Inside LEFT rarefaction fan
        // Riemann invariant R+ = u + 2a/(γ-1) = const = uL + 2aL/(γ-1)
        // Self-similar: ξ = u − a  (left-going characteristic)
        // Solving: u_fan = 2/(γ+1) * (aL + (γ-1)/2*uL + ξ)
        const double u_fan = 2.0 / (g + 1.0)
                             * (aL + (g - 1.0) / 2.0 * m_L.u + xi);
        const double a_fan = aL + (g - 1.0) / 2.0 * (m_L.u - u_fan);
        const double p_fan = m_L.p   * std::pow(a_fan / aL, 2.0 * g / (g - 1.0));
        const double r_fan = m_L.rho * std::pow(a_fan / aL, 2.0 / (g - 1.0));
        return {x, r_fan, u_fan, p_fan};
    }
    else if(xi <= m_uStar)
    {
        // Region 3* — left star state
        return {x, rhoSL, m_uStar, m_pStar};
    }
    else if(xi <= SR_left)
    {
        // Region 4* — right star state
        return {x, rhoSR, m_uStar, m_pStar};
    }
    else if(xi <= SR_right)
    {
        // Inside RIGHT rarefaction fan
        // Riemann invariant R- = u - 2a/(γ-1) = const = uR - 2aR/(γ-1)
        // Self-similar: ξ = u + a  (right-going characteristic)
        // Solving: u_fan = 2/(γ+1) * (ξ − aR + (γ-1)/2*uR)
        const double u_fan = 2.0 / (g + 1.0)
                             * (xi - aR + (g - 1.0) / 2.0 * m_R.u);
        const double a_fan = aR + (g - 1.0) / 2.0 * (u_fan - m_R.u);
        const double p_fan = m_R.p   * std::pow(a_fan / aR, 2.0 * g / (g - 1.0));
        const double r_fan = m_R.rho * std::pow(a_fan / aR, 2.0 / (g - 1.0));
        return {x, r_fan, u_fan, p_fan};
    }
    else
    {
        // Region 5 — undisturbed right state
        return {x, m_R.rho, m_R.u, m_R.p};
    }
}

// ===========================================================================
//  sampleLine  —  N points at cell-centres across [0, Lx]
// ===========================================================================

std::vector<SodExact::Point> SodExact::sampleLine(int N, double Lx,
                                                  double t, double x0) const
{
    std::vector<Point> result;
    result.reserve(N);
    for(int i = 0; i < N; ++i)
    {
        double x = Lx * (i + 0.5) / static_cast<double>(N);
        result.push_back(sample(x, t, x0));
    }
    return result;
}
