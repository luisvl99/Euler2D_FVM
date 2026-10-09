"""Helpers for scripted studies with the Euler2D_FVM command-line runner.

    io        read the CSV files and run.json written by euler_cli / the GUI
    exact     exact reference solutions (Riemann problem, smooth wave)
    runner    run euler_cli with cached results
    metrics   L1 errors and observed orders of convergence
    plotting  fixed colours per configuration, figure and table output
"""
