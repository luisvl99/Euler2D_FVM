"""Run euler_cli and cache its output under analysis/runs/."""

import os
import subprocess
from itertools import product
from pathlib import Path

from .io import read_run

ANALYSIS = Path(__file__).resolve().parents[1]
ROOT = ANALYSIS.parent
RUNS = ANALYSIS / "runs"

# The 8 solver configurations, as euler_cli options
CONFIGS = [dict(flux=f, recon=r, time=t)
           for f, r, t in product(("rusanov", "hllc"), ("pc", "muscl"), ("fe", "rk2"))]


def cli_path():
    """Locate euler_cli: $EULER_CLI, else the repository's build folder."""
    if "EULER_CLI" in os.environ:
        return Path(os.environ["EULER_CLI"])
    for folder in ("build", "build/Release"):
        for name in ("euler_cli", "euler_cli.exe"):
            path = ROOT / folder / name
            if path.exists():
                return path
    raise FileNotFoundError(
        "euler_cli not found. Build it with\n"
        "  cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build\n"
        "or set EULER_CLI to its path.")


def run(case, **options):
    """Run one case and return (output folder, run.json summary).

    Keyword options map to euler_cli options: nx=200 -> --nx 200,
    t_end=0.1 -> --t-end 0.1, subsonic=True -> --subsonic.  Results are
    cached: a run whose run.json already exists is not repeated.  Delete
    analysis/runs/ after changing the solver.  A blow-up is not an error
    here; check summary["status"].
    """
    name = "_".join([case] + [f"{key}{value}" for key, value in sorted(options.items())])
    out = RUNS / name
    if not (out / "run.json").exists():
        args = [str(cli_path()), "--case", case, "--out", str(out)]
        for key, value in sorted(options.items()):
            flag = "--" + key.replace("_", "-")
            args += [flag] if value is True else [flag, str(value)]
        result = subprocess.run(args, capture_output=True, text=True)
        if result.returncode not in (0, 2):      # 2 = blow-up, recorded in run.json
            raise RuntimeError(f"euler_cli failed: {result.stderr.strip()}")
    return out, read_run(out)
