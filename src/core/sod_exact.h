#pragma once
#include <cmath>
#include <vector>

// ---------------------------------------------------------------------------
// SodExact
//
//  Exact solution of the 1-D Riemann problem for the Euler equations,
//  implemented via Newton iteration (Toro, "Riemann Solvers and Numerical
//  Methods for Fluid Dynamics", Chapter 4).
//
//  The solution has up to 5 regions separated by 3 waves:
//
//    1 (left) | left wave | 3* | contact | 4* | right wave | 5 (right)
//
//  For the standard Sod problem the left wave is a rarefaction fan and
//  the right wave is a shock.  The implementation handles all combinations:
//  shock or rarefaction on each side.
//
//  Usage
//  -----
//    SodExact::State L{1.0, 0.0, 1.0};    // rho, u, p
//    SodExact::State R{0.125, 0.0, 0.1};
//    SodExact ex(L, R);
//    auto pts = ex.sampleLine(500, 1.0, 0.2, 0.5);  // N, Lx, t, x0
//
//  Expected Sod answer at t = 0.2:  p* ≈ 0.3031,  u* ≈ 0.9275
// ---------------------------------------------------------------------------

class SodExact
{
public:
    struct State { double rho = 0, u = 0, p = 0; };
    struct Point { double x = 0, rho = 0, u = 0, p = 0; };

    SodExact() = default;
    SodExact(State L, State R, double gamma = 1.4);

    // (Re)configure and solve — safe to call multiple times
    void setStates(State L, State R, double gamma = 1.4);

    // Sample exact solution at a single point (x, t),
    // initial discontinuity at x0
    Point sample(double x, double t, double x0) const;

    // Sample at N equally-spaced cell-centres over [0, Lx]
    std::vector<Point> sampleLine(int N, double Lx,
                                  double t, double x0) const;

    bool   isValid() const { return m_solved; }
    double pStar()   const { return m_pStar;  }
    double uStar()   const { return m_uStar;  }

private:
    State  m_L, m_R;
    double m_gamma  = 1.4;
    double m_pStar  = 0.0;
    double m_uStar  = 0.0;
    bool   m_solved = false;

    // Pressure function f(p, K) and its derivative — used in Newton iteration
    double f   (double p, const State& K, double aK) const;
    double dfdp(double p, const State& K, double aK) const;

    // Density in the star region adjacent to state K
    double starDensity(double pStar, const State& K) const;

    void solve();   // Newton loop → fills m_pStar, m_uStar

    static double soundSpeed(const State& K, double g)
    { return std::sqrt(g * K.p / K.rho); }
};
