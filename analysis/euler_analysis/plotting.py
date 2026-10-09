"""Shared figure and table output.

Every configuration keeps one colour in every figure.  Figures go to
docs/figures/, tables to docs/results/.
"""

import matplotlib

matplotlib.use("Agg")                     # write files, never open a window
import matplotlib.pyplot as plt           # noqa: E402

from .runner import ROOT                  # noqa: E402

FIGURES = ROOT / "docs" / "figures"
RESULTS = ROOT / "docs" / "results"

_NAMES = {"rusanov": "Rusanov", "hllc": "HLLC", "pc": "PC", "muscl": "MUSCL",
          "fe": "FE", "rk2": "SSP-RK2"}

_COLORS = {
    ("hllc", "pc", "fe"): "#e6194b",       # red
    ("hllc", "pc", "rk2"): "#4363d8",      # blue
    ("hllc", "muscl", "fe"): "#3cb44b",    # green
    ("hllc", "muscl", "rk2"): "#f58231",   # orange
    ("rusanov", "pc", "fe"): "#911eb4",    # purple
    ("rusanov", "pc", "rk2"): "#42d4f4",   # cyan
    ("rusanov", "muscl", "fe"): "#f032e6",  # magenta
    ("rusanov", "muscl", "rk2"): "#808000",  # olive
}

plt.rcParams.update({"axes.grid": True, "grid.alpha": 0.3, "font.size": 9,
                     "legend.fontsize": 8, "figure.dpi": 100})


def label(cfg):
    """'HLLC / MUSCL / SSP-RK2' for a configuration dict."""
    return " / ".join(_NAMES[cfg[k]] for k in ("flux", "recon", "time"))


def color(cfg):
    return _COLORS[(cfg["flux"], cfg["recon"], cfg["time"])]


def save_figure(fig, name):
    FIGURES.mkdir(parents=True, exist_ok=True)
    path = FIGURES / name
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {path.relative_to(ROOT)}")


def write_table(name, df, floatfmt="{:.6f}"):
    """Write df as docs/results/<name>.csv and as a Markdown table <name>.md."""
    RESULTS.mkdir(parents=True, exist_ok=True)
    df.to_csv(RESULTS / f"{name}.csv", index=False)

    def cell(value):
        return floatfmt.format(value) if isinstance(value, float) else str(value)

    lines = ["| " + " | ".join(df.columns) + " |",
             "|" + "|".join("---" for _ in df.columns) + "|"]
    lines += ["| " + " | ".join(cell(v) for v in row) + " |"
              for row in df.itertuples(index=False)]
    (RESULTS / f"{name}.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {(RESULTS / name).relative_to(ROOT)}.csv/.md")
    print("\n".join(lines))
