# Analysis scripts

Python scripts that run `euler_cli` over parameter grids and produce the figures in [`docs/figures`](../docs/figures) and the tables in [`docs/results`](../docs/results). Every number in the main README comes from here.

## Setup

Build `euler_cli` in Release mode (Debug builds are much slower). From the repository root:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DEULER_BUILD_GUI=OFF
```

```bash
cmake --build build
```

Then install the Python dependencies (Python 3.9 or newer):

```bash
pip install -r analysis/requirements.txt
```

The scripts look for `build/euler_cli` (or `build/Release/euler_cli` with multi-configuration generators). Set the environment variable `EULER_CLI` to use another executable.

## Running

Run every study (about four minutes from scratch):

```bash
python analysis/run_all.py
```

Or one study, for example:

```bash
python analysis/sod_schemes.py
```

Runs are cached in `analysis/runs/` (not versioned): a run whose `run.json` exists is not repeated. **Delete `analysis/runs/` after changing the solver**, otherwise the figures show old results.

## Studies

| Script | Question | Output |
|---|---|---|
| `sod_schemes.py` | Which of the 8 flux / reconstruction / time-scheme combinations is most accurate on Sod (N = 200, CFL 0.5)? | `sod_schemes.png`, `sod_schemes.md` |
| `sod_convergence.py` | How does the L1 error fall with N when the solution has discontinuities, and what does each accuracy level cost? | `sod_convergence.png`, `sod_cost.png`, `sod_convergence.md` |
| `cfl_sweep.py` | How do accuracy and robustness depend on the CFL number? | `cfl_sweep.png`, `cfl_sweep.md` |
| `smooth_order.py` | Do the schemes reach their design order on a smooth solution? | `smooth_order.png`, `smooth_order.md` |
| `shu_osher.py` | How well does each scheme keep the entropy waves behind a Mach-3 shock? | `shu_osher.png`, `shu_osher.md` |
| `riemann2d.py` | What do the 2-D Riemann problems look like with a diffusive and an accurate scheme? | `riemann2d.png` |
| `channel.py` | Does the channel flow reach a steady state with each boundary-condition set? | `channel.png`, `channel.md` |

## Package

`euler_analysis/` holds the shared pieces:

- `runner.run(case, **options)` runs `euler_cli` (options map to command-line flags, e.g. `nx=200`, `t_end=0.1`, `subsonic=True`) and returns the output folder and the `run.json` summary. `runner.CONFIGS` lists the 8 scheme combinations.
- `io.read_csv(path)` returns a pandas DataFrame and the `#` metadata header as a dict.
- `exact.riemann(x, t)` is the exact Riemann solution (Sod by default); `exact.smooth_wave(x, t)` is the smooth-wave solution. Run `python -m euler_analysis.exact` from this folder for a self-check.
- `metrics.l1_error` and `metrics.observed_orders`.
- `plotting` keeps one colour per configuration across all figures and writes figures and tables.

The files written by the GUI export have the same format, so `io.read_csv` also reads them.
