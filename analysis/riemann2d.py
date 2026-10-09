"""S6 - 2-D Riemann problems (Lax & Liu configurations 3, 4, 6 and 12).

256 x 256 cells on the unit square, CFL 0.5, each configuration's end time.
Compares the most diffusive scheme (Rusanov / PC) with the most accurate
one (HLLC / MUSCL), both with SSP-RK2, using the same 30 density contour
levels in each column.  Writes docs/figures/riemann2d.png.
"""

import numpy as np
import matplotlib.pyplot as plt

from euler_analysis import io, plotting, runner

N = 256
CASES = ["riemann3", "riemann4", "riemann6", "riemann12"]
SCHEMES = [dict(flux="rusanov", recon="pc", time="rk2"),
           dict(flux="hllc", recon="muscl", time="rk2")]


def density(case, cfg):
    folder, run = runner.run(case, nx=N, ny=N, **cfg)
    df, _ = io.read_csv(folder / "final.csv")
    shape = (N, N)                         # rows = j (y), columns = i (x)
    return (df.cx.to_numpy().reshape(shape), df.cy.to_numpy().reshape(shape),
            df.rho.to_numpy().reshape(shape), run["sim_time"])


def main():
    fig, axes = plt.subplots(2, len(CASES), figsize=(14, 7.4))
    for col, case in enumerate(CASES):
        fields = [density(case, cfg) for cfg in SCHEMES]
        rho_ref = fields[-1][2]
        levels = np.linspace(rho_ref.min(), rho_ref.max(), 30)
        for row, ((x, y, rho, t), cfg) in enumerate(zip(fields, SCHEMES)):
            ax = axes[row, col]
            ax.contour(x, y, rho, levels=levels, colors="k", linewidths=0.4)
            ax.set_aspect("equal")
            ax.set_xticks([0, 0.5, 1])
            ax.set_yticks([0, 0.5, 1])
            ax.grid(False)
            ax.set_title(f"config {case[7:]}, t = {t:g}\n{plotting.label(cfg)}", fontsize=8)
    fig.suptitle(f"2-D Riemann problems, density, {N} x {N} cells")
    plotting.save_figure(fig, "riemann2d.png")


if __name__ == "__main__":
    main()
