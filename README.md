# U(2) BPS counting and the nonlinear constraint

Reproducible numerical workspace for the calculations developed while studying the U(2) specialization of arXiv:2512.19946, with particular emphasis on Eq. (2.73), Appendix A, the positive BPS degeneracies, the single-saddle index envelope, and the nonlinear constraint.

The repository intentionally separates **exact coefficient generation**, **single-saddle formulas**, and **post-processing / asymptotic diagnostics**.  The main lesson from the numerical work is that the direct Appendix-A test is a **fixed microscopic charge Q** test: vary `2J`, compare the exact positive degeneracy `d[J,Q]` with the single-saddle count `d_sp[2J+Q]`, and study the crossing.

## Conventions

We specialize

```text
p = q = x = exp(2 pi i tau),
y = exp(2 pi i alpha),
```

and use

```text
Z_U(2)(x,y) = Tr_BPS x^(J1+J2+Q) y^Q.
```

On the equal-spin slice `J1=J2=J`,

```text
j = 2J + Q.
```

Throughout this repository, phrases such as **“through `q^N`”** or **“q-power exponent `j<=N`”** refer to the highest exponent retained in the Taylor expansion.  We do **not** mean a numerical inequality on the fugacity `q` itself.  For example, “through `q^200`” means that coefficients of `q^j` are retained for `j<=200`.

To expose the one-third grading,

```text
x = r^3,  y = s^3,
R = 3j,   S = 3Q.
```

The U(2) adjoint weights in the relative gauge fugacity are `{0,0,+1,-1}`.

## Exact oscillator product

The positive-counting oscillator families are

| family | monomial | multiplicity | k |
|---|---|---:|---:|
| scalar boson | `r^(2+3k) s^2` | `3(k+1)` | `k>=0` |
| chiral fermion | `r^(4+3k) s` | `3(k+1)` | `k>=0` |
| vector boson | `r^(6+3k)` | `k+1` | `k>=0` |
| vector fermion | `r^(3k) s^3` | `k+1` | `k>=1` |

Each family is multiplied over the four adjoint weights.  Bosons contribute `(1-X)^(-m)` and fermions `(1+X)^m`.

Rather than expand the full product directly, `src/u2bps/exact_series.py` uses

```text
D = r d/dr log Z,
R Z_R = sum_{k=1}^R D_k Z_{R-k}.
```

The U(2) Haar projection is

```text
c_(R,S) = [z^0] Z_(R,S) - [z^1] Z_(R,S).
```

The implementation uses Python arbitrary-size integers.  It reproduces exactly the stored coefficients through the term $q^{50}$ (equivalently, q-power exponent $j\le 50$) and is the transparent reference implementation for the coefficient calculation.

## Single saddle and nonlinear constraint

For N=2, M=1, k=-1,

```text
F_-1(omega) = -(16/27) (omega-pi i)^3 / omega^2,
S_-1(j,Q) = ext_omega [F_-1(omega) - omega j + pi i Q].
```

With `omega=i*pi*x`, the saddle equation is

```text
(j+16/27) x^3 - (16/9) x + 32/27 = 0.
```

The root with largest positive `Re S` is used and

```text
log d_sp(j) = Re S_-1(j).
```

The continuous N=2 nonlinear-constraint curve is parametrized by

```text
Q_ctr(a) = 8a/(1-a)^2,
j(a)     = 8a(1+a)^2/(1-a)^3,
2J_ctr   = j-Q_ctr.
```

At fixed microscopic `Q`, the periodic images are

```text
Q_ctr = Q + 2m.
```

Selecting the nearest image **after** measuring the numerical crossing is a consistency diagnostic, not a prediction of the image label `m`.

## Reproduced numerical milestones

### Exact positive coefficients

- Through the coefficient of `q^50` (q-power exponent `j<=50`): 3,874 nonzero `(j,Q)` coefficients.
- At q-power exponent `j=50` (the coefficient of `q^50`), the fixed-`j` maximum is at `Q=17` with

```text
1499914493457031358970381412493.
```

- Through the coefficient of `q^100` (q-power exponent `j<=100`), the peak locations include

```text
j = 50,60,70,80,90,100
Q_peak = 17,20,22,24,26,28.
```

A fit over q-power exponents `j=50..100` gives approximately

```text
Q_peak ~ 1.291882 j^(2/3),
log c_max ~ 6.112350 j^(2/3) - 2.135663 log j - 5.118475.
```

This fixed-j peak analysis is useful characterization, but **is not the direct Appendix-A test**.

### Raw index control experiment

At fixed j, setting y=-1 gives an alternating sum over charge.  The exact raw index undergoes enormous cancellations.  This is a negative/control result: the raw `|index|` is not the same object as the one-saddle envelope `d_sp` used in Appendix A.

### Direct Appendix-A crossing

For fixed `Q`, solve

```text
d[J,Q] = d_sp[2J+Q].
```

The benchmark sector `Q=4/3` gives

```text
2J_cross = 22.003694,
Appendix-A quoted prediction = 21.9575,
relative difference ~ 0.21%.
```

For integer Q at larger spin, representative results include

```text
Q=3:  2J_cross=103.251119, nearest image 2J_ctr=103.776695
Q=6:  2J_cross=328.862367, nearest image 2J_ctr=328.302031
```

### Full one-third charge lattice

The calculation was extended to all

```text
Q in (1/3) Z
```

rather than only integer Q.  The three residue classes show similar behavior.  The crossing stays close, in a **relative** sense, to the family of periodic nonlinear-constraint images.

### High-Q oscillations

The reorganized variables

```text
n = 2J = (R-S)/3,
S = 3Q
```

allow a rectangular high-spin / low-charge computation.  Production runs reached `Q=40` and `2J` of order 6000 in the crossing data.

A crucial observation is that two notions of convergence behave differently:

```text
absolute residual: |2J_cross - 2J_ctr|  does not decrease,
relative residual: |2J_cross/2J_ctr - 1| does decrease.
```

Successive large-Q absolute residual maxima grow, while their relative size falls below the percent level.

### Ratio convergence

The cleanest macroscopic diagnostic is

```text
R_J(Q) = 2J_cross / 2J_ctr^closest.
```

Through Q=40 the ratio oscillates about 1 with decreasing relative amplitude.  For example

```text
Q=20:  R_J=1.009359
Q=30:  R_J=1.001286
Q=35:  R_J=0.996359
Q=40:  R_J=1.002469.
```

The supported asymptotic statement is therefore

```text
2J_cross / 2J_ctr^closest -> 1
```

with oscillatory corrections.  The data do **not** support the stronger statement `2J_cross-2J_ctr -> 0`.

### Vertical mismatch and shrinking slope

Define

```text
G(2J,Q) = log d[J,Q] - log d_sp[2J+Q].
```

At the analytic nonlinear-constraint point,

```text
Delta_log(Q) = G(2J_ctr,Q),
s(Q) = dG/d(2J)|_(2J_ctr).
```

The measured crossing shift satisfies extremely accurately

```text
delta(2J) = 2J_cross-2J_ctr ~ -Delta_log/s.
```

The slope scales numerically as

```text
|s| ~ 0.449 Q_ctr^(-0.495),
```

very close to `Q_ctr^(-1/2)`.  The vertical mismatch itself remains of order one over the directly measured range.  Therefore an O(1), or slowly varying/logarithmic, correction to `log d` can naturally generate an apparently power-like horizontal shift

```text
delta(2J) ~ sqrt(Q_ctr).
```

This is why the growing absolute horizontal displacement should not be interpreted directly as a power-law correction to the entropy.

### Gaussian prefactor diagnostic

The simple one-dimensional Gaussian phase of the final omega saddle,

```text
-1/2 arg S''(omega_*),
```

was tested separately.  It produces only a roughly 0.09--0.12 shift in charge units and does not remove the observed oscillations.  This should **not** be confused with the full matrix-integral one-loop determinant.

## Repository layout

```text
src/u2bps/
  exact_series.py   exact Python-integer Taylor recurrence + Haar projection
  saddle.py         single saddle, nonlinear constraint, periodic images
  crossings.py      fixed-Q crossings and vertical-mismatch diagnostics
  rectangular.py    readable high-spin / low-charge reference recurrence

scripts/
  01_generate_exact_series.py
  02_peak_locus.py
  03_raw_index_control.py
  04_appendixA_direct.py
  05_periodic_constraint_scan.py
  06_generate_rectangular_strip.py
  07_oscillation_scaling.py
  08_gaussian_prefactor_test.py
  09_vertical_mismatch.py
  10_ratio_convergence.py
  reproduce_analysis.py

data/reference/
  stored coefficient/crossing tables from the production runs

docs/
  detailed interpretation and chronology of the tests
```

## Installation

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e .
```

## Quick reproduction

The lightweight analyses use the stored production tables:

```bash
python scripts/reproduce_analysis.py
```

The exact coefficients through `q^50` (maximum q-power exponent `j=50`) can be regenerated from scratch with

```bash
python scripts/01_generate_exact_series.py --qmax 50 --progress
```

For a modest direct `(2J,3Q)` reference run:

```bash
python scripts/06_generate_rectangular_strip.py --nmax 200 --Qmax 6 --n-theta 64
```

Increase `n_theta` and compare results to test gauge-Fourier convergence.

## Important interpretation warnings

1. The fixed-j coefficient maximum `Q_peak(j)` is not the Appendix-A curve.
2. The raw exact index magnitude is not the one-saddle envelope because of large cancellations.
3. The nearest periodic image is an a-posteriori comparison unless a separate argument predicts the image label.
4. The simple final-omega Gaussian prefactor is not the full one-loop determinant.
5. A power-like horizontal crossing shift does not imply a power-like correction to `log d`; the crossing slope becomes small at large charge.

## Paper

The project was developed around arXiv:2512.19946.  Please cite the original paper when using the physical formulas; this repository is a numerical/reproducibility workspace rather than an independent derivation of the underlying index technology.
