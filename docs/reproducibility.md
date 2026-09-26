# Reproducibility tiers

The repository contains two kinds of computations.

## Tier 1: exact / lightweight

These can be rerun directly from the Python source on an ordinary workstation:

- exact sparse coefficient recurrence at modest q-power exponent `j`,
- single-saddle evaluation,
- nonlinear-constraint curve and periodic images,
- fixed-Q crossing extraction,
- raw-index control analysis,
- oscillation/ratio plots from stored crossing tables,
- Gaussian final-omega prefactor diagnostic,
- vertical mismatch and slope analysis.

The exact recurrence has been regression-tested against the stored coefficients through `q^50` (q-power exponent `j<=50`).

## Tier 2: production high-spin strips

The largest runs (Q<=40, 2J of order 6000) use the same `(2J,3Q)` reorganization described in `src/u2bps/rectangular.py`, together with banding/chunking and increased gauge-Fourier resolution.  Those production runs were computationally heavier, so the compact crossing tables are shipped under `data/reference/`.

The readable `rectangular.py` implementation is mathematically the same recurrence but is optimized for transparency rather than maximum throughput.  It has been checked against the exact sparse coefficients at low n and S.  For new high-spin production runs, increase `n_theta` and compare the resulting crossings before trusting the last digits.

## Regression checks performed while assembling this repository

- exact sparse recurrence versus stored coefficients through `q^50`: zero mismatches in the tested low-`j` window;
- rectangular `(2J,3Q)` Fourier implementation versus exact coefficients for `2J<=10`, `Q<=2`: zero mismatches at 64 gauge angles to numerical precision;
- Q=4/3 direct crossing: `2J=22.0036943161`;
- vertical-mismatch slope exponent: `-0.4954`.
