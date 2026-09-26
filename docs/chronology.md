# Chronology of the numerical investigation

This file records what was tested, why some apparently natural tests were abandoned, and which diagnostics should be regarded as physically meaningful.

## 1. Exact positive coefficient distribution

The first calculation expanded the positive U(2) BPS product and Haar-projected it exactly.  The fixed-`j` charge profile has a peak near `Q ~ j^(2/3)`, where `j` is the exponent of `q` in the Taylor expansion.  This was initially compared with the principal nonlinear-constraint curve.

That comparison is not the direct Appendix-A test because Appendix A fixes microscopic Q and scans angular momentum.

## 2. Raw y=-1 index

The exact alternating sum was computed at fixed j.  It is many orders of magnitude smaller than the positive degeneracy because of cancellations.  No exact discrete equality with a positive coefficient should be expected.  The paper's Appendix-A comparison removes index oscillations by isolating one saddle, so the raw index was retained only as a control experiment.

## 3. Single-saddle fixed-Q crossing

The correct direct observable is

`f_Q(2J)=log d_sp(2J+Q)/log d[J,Q]`.

A crossing `f_Q=1` defines `2J_cross`.  The Q=4/3 benchmark reproduces the quoted Appendix-A value at the 0.2% level.

## 4. Periodic images

The nonlinear constraint is modulo 2 in charge.  At fixed microscopic Q, one compares with representatives `Q_ctr=Q+2m`.  The nearest representative tracks the numerical crossing well at moderate and large charge, but choosing m after seeing the data is a consistency check, not a branch-selection theory.

## 5. Full Q in one-third units

The integer-only high-spin calculation was generalized to every S=3Q.  All three residue classes behave similarly.  The direct wrapped mod-2 residual continues to oscillate by O(1); there is no simple monotonic phase convergence.

## 6. Higher Q and J

Reorganizing the grading in terms of `(2J,3Q)` made large-J fixed-Q calculations much cheaper.  Crossings were followed through Q=40, reaching 2J of order 6000.

The absolute residual grows at some successive peaks, while the relative residual shrinks.  This motivates using `2J_cross/2J_ctr` as the clean macroscopic convergence diagnostic.

## 7. One-dimensional Gaussian prefactor

The phase of the final omega Gaussian was added.  It is much too small to explain the order-one wrapped phase excursions and does not systematically improve the crossing residual.  This is not a calculation of the full matrix-integral one-loop determinant.

## 8. Vertical mismatch

The apparent `sqrt(Q)` horizontal residual was traced to the fact that the crossing function gets flatter:

`|d/d(2J)(log d-log d_sp)| ~ Q_ctr^(-1/2)`.

The vertical mismatch in log degeneracy is O(1) over the directly resolved range, and

`delta(2J) ~ -Delta_log/slope`

reconstructs the measured crossing almost exactly.  Therefore the horizontal power law is primarily an amplification effect and should not be called a power-law entropy correction.
