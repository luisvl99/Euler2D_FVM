"""S5 - Shu-Osher: how well each configuration keeps the entropy waves behind
the shock.

N = 400 on the canonical domain (Lx = 10, plotted as x in [-5, 5]), CFL 0.5,
t = 1.8.  There is no exact solution; the reference is HLLC / MUSCL /
SSP-RK2 on N = 6400, averaged onto the coarse cells for the L1 error.
Writes docs/figures/shu_osher.png and docs/results/shu_osher.{csv,md}.
"""

import pandas as pd
import matplotlib.pyplot as plt

from euler_analysis import io, metrics, plotting, runner

N, N_REF = 400, 6400


def main():
    folder, _ = runner.run("shu-osher", nx=N_REF, flux="hllc", recon="muscl", time="rk2")
    ref, _ = io.read_csv(folder / "probe_j00.csv")
    ref_on_coarse = ref.rho.to_numpy().reshape(N, -1).mean(axis=1)

    fig, (full, zoom) = plt.subplots(1, 2, figsize=(11, 4.2))
    for ax in (full, zoom):
        ax.plot(ref.x - 5.0, ref.rho, "k-", lw=1.2, label=f"Reference (N = {N_REF})")

    rows = []
    for cfg in runner.CONFIGS:
        folder, _ = runner.run("shu-osher", nx=N, **cfg)
        df, _ = io.read_csv(folder / "probe_j00.csv")
        rows.append({"Configuration": plotting.label(cfg),
                     "L1(rho) vs reference": metrics.l1_error(df.rho, ref_on_coarse)})
        for ax in (full, zoom):
            ax.plot(df.x - 5.0, df.rho, "-", lw=0.9, color=plotting.color(cfg),
                    label=plotting.label(cfg))

    full.set(xlabel="x", ylabel="density", title=f"Shu-Osher, t = 1.8, N = {N}")
    zoom.set(xlabel="x", xlim=(0.4, 2.6), ylim=(2.8, 4.9),
             title="Entropy waves behind the shock")
    full.legend(loc="lower left")
    plotting.save_figure(fig, "shu_osher.png")

    table = pd.DataFrame(rows).sort_values("L1(rho) vs reference", ascending=False)
    plotting.write_table("shu_osher", table, floatfmt="{:.4f}")


if __name__ == "__main__":
    main()
