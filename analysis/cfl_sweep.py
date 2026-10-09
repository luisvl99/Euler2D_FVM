"""S3 - Sod shock tube: accuracy and robustness against the CFL number.

N = 200, t = 0.2, CFL = 0.1 ... 3.0.  A run "fails" when the solver stops
with a non-physical state (exit code 2).  Remember that this CFL is about
twice the classical Courant number (see the README).
Writes docs/figures/cfl_sweep.png and docs/results/cfl_sweep.{csv,md}.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from euler_analysis import exact, io, metrics, plotting, runner

N = 200
CFLS = np.round(np.arange(0.1, 3.0001, 0.1), 2)


def main():
    fig, ax = plt.subplots(figsize=(7, 4.5))
    rows = []

    for cfg in runner.CONFIGS:
        errors, last_ok = [], None
        for cfl in CFLS:
            folder, run = runner.run("sod", nx=N, cfl=cfl, **cfg)
            if run["status"] == "blowup":
                errors.append(np.nan)
                continue
            if all(np.isfinite(errors)):     # no failure at a smaller CFL yet
                last_ok = cfl
            df, _ = io.read_csv(folder / "probe_j00.csv")
            errors.append(metrics.l1_error(df.rho, exact.riemann(df.x, run["sim_time"])[0]))

        best = int(np.nanargmin(errors))
        rows.append({"Configuration": plotting.label(cfg),
                     "largest CFL without failure": last_ok,
                     "most accurate CFL": CFLS[best],
                     "L1(rho) there": errors[best]})
        ax.semilogy(CFLS, errors, "o-", ms=3, color=plotting.color(cfg), label=plotting.label(cfg))

    ax.axvline(2.0, color="k", ls=":", lw=0.8)
    ax.text(2.02, ax.get_ylim()[1] * 0.8, "classical\nCourant = 1", fontsize=7, va="top")
    ax.set(xlabel="CFL (GUI / solver definition)", ylabel="L1 error of density",
           title=f"Sod: accuracy against CFL (N = {N}); gaps = solver stopped")
    ax.legend(ncol=2)
    plotting.save_figure(fig, "cfl_sweep.png")

    plotting.write_table("cfl_sweep", pd.DataFrame(rows), floatfmt="{:.4g}")


if __name__ == "__main__":
    main()
