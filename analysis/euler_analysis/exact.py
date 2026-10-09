"""Exact reference solutions.

riemann()      1-D Riemann problem of the Euler equations (Toro, ch. 4), a
               port of the C++ SodExact class; defaults to Sod's shock tube.
smooth_wave()  density of the smooth-wave test case.
"""

import numpy as np

GAMMA = 1.4
SOD_LEFT = (1.0, 0.0, 1.0)      # (rho, u, p)
SOD_RIGHT = (0.125, 0.0, 0.1)


def _pressure_function(p, rho, pk, a, g):
    """f_K(p) and its derivative (Toro eq. 4.6 and 4.37)."""
    if p > pk:                                   # shock
        A = 2.0 / ((g + 1.0) * rho)
        B = (g - 1.0) / (g + 1.0) * pk
        q = np.sqrt(A / (p + B))
        return (p - pk) * q, q * (1.0 - (p - pk) / (2.0 * (p + B)))
    ratio = p / pk                               # rarefaction
    return (2.0 * a / (g - 1.0) * (ratio ** ((g - 1.0) / (2.0 * g)) - 1.0),
            ratio ** (-(g + 1.0) / (2.0 * g)) / (rho * a))


def star_state(left=SOD_LEFT, right=SOD_RIGHT, gamma=GAMMA, tol=1e-12, max_iter=100):
    """Pressure and velocity (p*, u*) of the star region, by Newton iteration."""
    (rl, ul, pl), (rr, ur, pr) = left, right
    al, ar = np.sqrt(gamma * pl / rl), np.sqrt(gamma * pr / rr)

    # PVRS initial guess (Toro eq. 9.28)
    p = max(0.5 * (pl + pr) - 0.125 * (ur - ul) * (rl + rr) * (al + ar), 1e-10)
    for _ in range(max_iter):
        fl, dl = _pressure_function(p, rl, pl, al, gamma)
        fr, dr = _pressure_function(p, rr, pr, ar, gamma)
        p_new = max(p - (fl + fr + ur - ul) / (dl + dr), 1e-10)
        converged = abs(p_new - p) < tol * 0.5 * (p_new + p)
        p = p_new
        if converged:
            break
    else:
        raise RuntimeError("exact Riemann solver: Newton iteration did not converge")

    fl, _ = _pressure_function(p, rl, pl, al, gamma)
    fr, _ = _pressure_function(p, rr, pr, ar, gamma)
    return p, 0.5 * (ul + ur) + 0.5 * (fr - fl)


def riemann(x, t, x0=0.5, left=SOD_LEFT, right=SOD_RIGHT, gamma=GAMMA):
    """Exact (rho, u, p) at positions x and time t > 0, discontinuity at x0."""
    g = gamma
    x = np.asarray(x, dtype=float)
    xi = (x - x0) / t
    ps, us = star_state(left, right, g)
    rho, u, p = np.empty_like(x), np.empty_like(x), np.empty_like(x)

    # s = -1 for the left wave, +1 for the right wave (Toro's symmetric form)
    for mask, (rk, uk, pk), s in ((xi <= us, left, -1.0), (xi > us, right, 1.0)):
        ak = np.sqrt(g * pk / rk)
        xs = xi[mask]
        if ps > pk:                                                    # shock
            speed = uk + s * ak * np.sqrt((g + 1) / (2 * g) * ps / pk + (g - 1) / (2 * g))
            mu = (g - 1) / (g + 1)
            outside = s * xs > s * speed
            rho[mask] = np.where(outside, rk, rk * (ps / pk + mu) / (mu * ps / pk + 1))
            u[mask] = np.where(outside, uk, us)
            p[mask] = np.where(outside, pk, ps)
        else:                                                          # rarefaction
            a_star = ak * (ps / pk) ** ((g - 1) / (2 * g))
            outside = s * xs >= s * (uk + s * ak)                      # beyond the head
            star = s * xs <= s * (us + s * a_star)                     # behind the tail
            u_fan = 2 / (g + 1) * (-s * ak + (g - 1) / 2 * uk + xs)
            a_fan = ak + s * (g - 1) / 2 * (u_fan - uk)
            with np.errstate(invalid="ignore"):                        # a_fan < 0 off the fan
                rho[mask] = np.select([outside, star], [rk, rk * (ps / pk) ** (1 / g)],
                                      rk * (a_fan / ak) ** (2 / (g - 1)))
                p[mask] = np.select([outside, star], [pk, ps],
                                    pk * (a_fan / ak) ** (2 * g / (g - 1)))
            u[mask] = np.select([outside, star], [uk, us], u_fan)
    return rho, u, p


def smooth_wave(x, t, lx=1.0):
    """Exact density of the smooth-wave case (see physics_testcases.h)."""
    s = (np.asarray(x, dtype=float) - 2.0 * t - 0.3 * lx) / (0.08 * lx)
    return 1.0 + 0.2 * np.exp(-s * s)


if __name__ == "__main__":
    # Self-check against the textbook Sod values (Toro) and mirror symmetry
    p_star, u_star = star_state()
    assert abs(p_star - 0.30313) < 1e-4 and abs(u_star - 0.92745) < 1e-4, (p_star, u_star)
    x = np.linspace(0.0, 1.0, 101)
    a = riemann(x, 0.2)
    b = riemann(1.0 - x, 0.2, left=SOD_RIGHT, right=SOD_LEFT)
    assert np.allclose(a[0], b[0]) and np.allclose(a[1], -b[1]) and np.allclose(a[2], b[2])
    print(f"exact.py OK: p* = {p_star:.5f}, u* = {u_star:.5f}")
