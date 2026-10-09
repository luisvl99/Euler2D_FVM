"""S2 + S8 - Sod shock tube: mesh convergence and cost of the 8 configurations.

N = 50 ... 1600 cells on Lx = 1, CFL = 0.5, t = 0.2.  The wall time comes
from run.json (solver loop only); use a Release build of euler_cli.
Writes docs/figures/sod_convergence.png, docs/figures/sod_cost.png and
docs/results/sod_convergence.{csv,md}.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from euler_analysis import exact, io, metrics, plotting, runner

SIZES = [50, 100, 200, 400, 800, 1600]


def main():
    fig, ax = plt.subplots(figsize=(8, 4.5))
    cost_fig, cost_ax = plt.subplots(figsize=(6, 4.5))
    rows = []

    for cfg in runner.CONFIGS:
        errors, walls = [], []
        for n in SIZES:
            folder, run = runner.run("sod", nx=n, cfl=0.5, **cfg)
            df, _ = io.read_csv(folder / "probe_j00.csv")
            errors.append(metrics.l1_error(df.rho, exact.riemann(df.x, run["sim_time"])[0]))
            walls.append(run["wall_time_s"])

        h = 1.0 / np.array(SIZES)
        orders = metrics.observed_orders(h, errors)
        rows.append({"Configuration": plotting.label(cfg),
                     **{f"N={n}": e for n, e in zip(SIZES, errors)},
                     "order (last 3 pairs)": float(np.mean(orders[-3:]))})
        style = dict(color=plotting.color(cfg), label=plotting.label(cfg), marker="o", ms=3)
        ax.loglog(SIZES, errors, **style)
        cost_ax.loglog(walls, errors, **style)

    # Reference slopes, anchored at the coarsest mesh
    n = np.array(SIZES, dtype=float)
    for slope, ls in ((0.5, ":"), (1.0, "--")):
        ax.loglog(n, 0.04 * (n / n[0]) ** -slope, "k" + ls, lw=0.8, label=f"slope {slope:g}")

    ax.set(xlabel="cells N", ylabel="L1 error of density",
           title="Sod: mesh convergence (CFL 0.5, t = 0.2)")
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1))
    plotting.save_figure(fig, "sod_convergence.png")

    cost_ax.set(xlabel="wall time of the solver loop [s]", ylabel="L1 error of density",
                title=f"Sod: accuracy against cost (N = {SIZES[0]} ... {SIZES[-1]})")
    cost_ax.legend()
    plotting.save_figure(cost_fig, "sod_cost.png")

    plotting.write_table("sod_convergence", pd.DataFrame(rows), floatfmt="{:.3g}")


if __name__ == "__main__":
    main()
