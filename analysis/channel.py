"""S7 - Channel flow: convergence to a steady state with each boundary set.

100 x 20 cells (1 m x 0.2 m), HLLC / MUSCL / SSP-RK2, CFL 0.5, starting at
rest with rho = 1.225 kg/m^3 and p = 101325 Pa.
Writes docs/figures/channel.png and docs/results/channel.{csv,md}.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from euler_analysis import io, plotting, runner

STEPS = 20000
RUNS = [
    ("supersonic, u = 600 m/s", dict(u=600)),
    ("subsonic, u = 68 m/s, p_back = p", dict(u=68, subsonic=True)),
    ("subsonic, u = 68 m/s, p_back = 0.99 p", dict(u=68, subsonic=True, p_back=100312)),
]


def column_mean(df, i, values):
    return float(values[df.i == i].mean())


def main():
    fig, ax = plt.subplots(figsize=(7, 4.2))
    rows = []
    for name, options in RUNS:
        folder, run = runner.run("channel", max_iter=STEPS, **options)
        res, _ = io.read_csv(folder / "residuals.csv")
        ax.semilogy(res.iter, res.residual_L2 / res.residual_L2.iloc[0], lw=1, label=name)

        df, meta = io.read_csv(folder / "final.csv")
        last = meta["Nx"] - 1
        mass = df.rho * df.u
        rows.append({"Run": name,
                     "residual drop (decades)":
                         float(np.log10(res.residual_L2.iloc[0] / res.residual_L2.iloc[-1])),
                     "u inlet": column_mean(df, 0, df.u), "u outlet": column_mean(df, last, df.u),
                     "p outlet": column_mean(df, last, df.p),
                     "mass flux in": column_mean(df, 0, mass),
                     "mass flux out": column_mean(df, last, mass)})

    ax.set(xlabel="iteration", ylabel="||R|| / ||R||(0)",
           title="Channel flow from rest (100 x 20, HLLC / MUSCL / SSP-RK2)")
    ax.legend()
    plotting.save_figure(fig, "channel.png")
    plotting.write_table("channel", pd.DataFrame(rows), floatfmt="{:.2f}")


if __name__ == "__main__":
    main()
