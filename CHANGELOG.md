# Changelog

## 1.0.0 (2026-10-09)

First public release.

- **Solver** (`src/core`): finite-volume method on a structured Cartesian mesh with Rusanov and HLLC fluxes, piecewise-constant and MUSCL (MinMod) reconstruction, Forward Euler and SSP-RK2 time integration, and ghost-cell boundary conditions (slip wall, symmetry, supersonic and characteristic subsonic inlet/outlet, far-field, zero gradient). The run stops with a clear message when the solution becomes non-physical.
- **Test cases**: channel flow, Sod shock tube, Shu–Osher, 2-D Riemann problems (Lax & Liu configurations 3, 4, 6, 12) and a smooth wave for order-of-accuracy studies.
- **GUI** (Qt 6): colour maps with zoom and pan, live residual plot, live 1-D profile viewer with the exact solution for Sod, CSV export with a metadata header.
- **Command-line runner** `euler_cli` for scripted studies.
- **Analysis scripts** (`analysis/`) that reproduce every figure and table in the README.
- **Tests**: headless regression checks (exact solver, flux properties, conservation, symmetry, order of accuracy, robustness) and CLI smoke tests, run by CI on Linux, Windows and macOS.
- **Windows binaries**: the GUI and `euler_cli` with the Qt runtime, attached to the GitHub release.
