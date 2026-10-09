"""S1 - Sod shock tube: the 8 solver configurations against the exact solution.

N = 200 cells on Lx = 1, CFL = 0.5, t = 0.2.
Writes docs/figures/sod_schemes.png and docs/results/sod_schemes.{csv,md}.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from euler_analysis import exact, io, metrics, plotting, runner

N, CFL = 200, 0.5


def main():
    fig, (full, zoom) = plt.subplots(1, 2, figsize=(11, 4.2))
    x = np.linspace(0.0, 1.0, 2000)
    for ax in (full, zoom):
        ax.plot(x, exact.riemann(x, 0.2)[0], "k-", lw=1.5, label="Exact")

    rows = []
    for cfg in runner.CONFIGS:
        folder, run = runner.run("sod", nx=N, cfl=CFL, **cfg)
        df, _ = io.read_csv(folder / "probe_j00.csv")
        rho_e, _, p_e = exact.riemann(df.x, run["sim_time"])
        rows.append({"Flux": plotting.label(cfg).split(" / ")[0],
                     "Reconstruction": cfg["recon"].upper(),
                     "Time scheme": "SSP-RK2" if cfg["time"] == "rk2" else "Forward Euler",
                     "L1(rho)": metrics.l1_error(df.rho, rho_e),
                     "L1(p)": metrics.l1_error(df.p, p_e)})
        for ax in (full, zoom):
            ax.plot(df.x, df.rho, "-", lw=1, color=plotting.color(cfg), label=plotting.label(cfg))

    full.set(xlabel="x", ylabel="density", title=f"Sod shock tube, t = 0.2, N = {N}, CFL = {CFL}")
    zoom.set(xlabel="x", xlim=(0.62, 0.9), ylim=(0.1, 0.45), title="Contact discontinuity and shock")
    zoom.legend(loc="upper right")
    plotting.save_figure(fig, "sod_schemes.png")

    table = pd.DataFrame(rows).sort_values("L1(rho)", ascending=False)
    plotting.write_table("sod_schemes", table)


if __name__ == "__main__":
    main()
