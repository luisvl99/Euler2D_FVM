"""Run every study and regenerate all figures and tables (a few minutes with a
Release build; cached runs in analysis/runs/ are reused)."""

import channel
import cfl_sweep
import riemann2d
import shu_osher
import smooth_order
import sod_convergence
import sod_schemes

for study in (sod_schemes, sod_convergence, cfl_sweep, smooth_order, shu_osher, riemann2d,
              channel):
    print(f"\n=== {study.__name__} ===")
    study.main()
