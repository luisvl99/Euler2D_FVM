"""S4 - Order of accuracy on a smooth solution (smooth-wave case).

A Gaussian density bump advected at u = 2; exact solution
exact.smooth_wave().  N = 50 ... 1600, CFL = 0.5, t = 0.2.
Writes docs/figures/smooth_order.png and docs/results/smooth_order.{csv,md}.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from euler_analysis import exact, io, metrics, plotting, runner

SIZES = [50, 100, 200, 400, 800, 1600]


def main():
    fig, ax = plt.subplots(figsize=(8, 4.5))
    rows = []

    for cfg in runner.CONFIGS:
        errors = []
        for n in SIZES:
            folder, run = runner.run("smooth-wave", nx=n, cfl=0.5, **cfg)
            df, _ = io.read_csv(folder / "probe_j00.csv")
            errors.append(metrics.l1_error(df.rho, exact.smooth_wave(df.x, run["sim_time"])))

        orders = metrics.observed_orders(1.0 / np.array(SIZES), errors)
        rows.append({"Configuration": plotting.label(cfg),
                     **{f"N={n}": e for n, e in zip(SIZES, errors)},
                     "order (800 -> 1600)": float(orders[-1])})
        ax.loglog(SIZES, errors, "o-", ms=3, color=plotting.color(cfg), label=plotting.label(cfg))

    n = np.array(SIZES, dtype=float)
    for slope, ls in ((1.0, "--"), (2.0, "-.")):
        ax.loglog(n, 0.02 * (n / n[0]) ** -slope, "k" + ls, lw=0.8, label=f"slope {slope:g}")

    ax.set(xlabel="cells N", ylabel="L1 error of density",
           title="Smooth wave: order of accuracy (CFL 0.5)")
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1))
    plotting.save_figure(fig, "smooth_order.png")

    plotting.write_table("smooth_order", pd.DataFrame(rows), floatfmt="{:.3g}")


if __name__ == "__main__":
    main()
