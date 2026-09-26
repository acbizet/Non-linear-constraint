#!/usr/bin/env python3
"""Run the lightweight analysis/plot scripts from the shipped reference data."""

import subprocess, sys
scripts=[
 "02_peak_locus.py",
 "03_raw_index_control.py",
 "04_appendixA_direct.py",
 "05_periodic_constraint_scan.py",
 "07_oscillation_scaling.py",
 "08_gaussian_prefactor_test.py",
 "09_vertical_mismatch.py",
 "10_ratio_convergence.py",
]
for s in scripts:
    print("\n===",s,"===")
    subprocess.run([sys.executable,"scripts/"+s],check=True)
