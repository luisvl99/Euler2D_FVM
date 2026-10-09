"""Error norms and convergence rates."""

import numpy as np


def l1_error(numeric, exact):
    """Mean absolute error over the cells.

    On a uniform mesh this is the integral of |error| divided by the domain
    length, so it does not depend on the number of cells.
    """
    return float(np.mean(np.abs(np.asarray(numeric) - np.asarray(exact))))


def observed_orders(h, errors):
    """Rates log(e_k / e_k+1) / log(h_k / h_k+1) between successive meshes."""
    h = np.asarray(h, dtype=float)
    e = np.asarray(errors, dtype=float)
    return np.log(e[:-1] / e[1:]) / np.log(h[:-1] / h[1:])
