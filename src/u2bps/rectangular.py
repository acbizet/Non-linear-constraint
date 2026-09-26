"""Reference high-spin / low-charge implementation in (2J,3Q) variables.

The exact sparse z-Laurent recurrence is ideal when the q-power exponent
j is at most O(100), but it wastes work
when the question is fixed, relatively small Q and very large J.  In that regime
use

    n = 2J = (R-S)/3,   S = 3Q.

The two zero-n oscillator families are absorbed into Z_0(S,z).  The remaining
product can again be generated from a logarithmic-derivative recurrence in n.

For very large n the gauge projection is performed by trapezoidal Fourier
quadrature over the relative U(2) angle.  Increasing `n_theta` gives a direct
convergence check.  This module is intentionally a readable reference
implementation.  The Q<=40 production runs used banded/chunked variants of the
same recurrence and checked stability under increased Fourier resolution.

Because a complete n~6000, S~120 run is CPU intensive, the repository ships the
small summary crossing tables and recommends first validating at n<=300.
"""

from __future__ import annotations

import math
import numpy as np
from numba import njit


def _poly_mul_truncated(a, terms, Smax):
    out = np.zeros_like(a)
    for s0 in range(Smax + 1):
        if a[s0] == 0:
            continue
        for ds, c in terms:
            if s0 + ds <= Smax:
                out[s0 + ds] += a[s0] * c
    return out


def zero_spin_base(Smax: int, theta: float) -> np.ndarray:
    """Z_{n=0}(S,theta) from scalar k=0 and vector-fermion k=1 modes."""
    z = np.exp(1j * theta)
    a = np.zeros(Smax + 1, dtype=np.complex128)
    a[0] = 1.0
    for w in (0, 0, 1, -1):
        zw = z**w
        # (1-u^2 z^w)^(-3)
        terms = []
        for ell in range(Smax // 2 + 1):
            terms.append((2 * ell, math.comb(ell + 2, 2) * zw**ell))
        a = _poly_mul_truncated(a, terms, Smax)
        # (1+u^3 z^w)^2
        terms = [(0, 1.0), (3, 2.0 * zw), (6, zw**2)]
        a = _poly_mul_truncated(a, terms, Smax)
    return a


def build_D_theta(nmax: int, Smax: int, theta: float) -> np.ndarray:
    """D[n,S] = coefficient of t^n u^S in t d/dt log Z at fixed theta."""
    D = np.zeros((nmax + 1, Smax + 1), dtype=np.float64)

    def add_mode(p, b, mult, fermion):
        if p <= 0:
            return
        hmax = nmax // p
        if b:
            hmax = min(hmax, Smax // b)
        for h in range(1, hmax + 1):
            sign = -1.0 if (fermion and h % 2 == 0) else 1.0
            adj = 2.0 + 2.0 * math.cos(h * theta)
            D[p * h, b * h] += mult * p * sign * adj

    # scalar boson: p=k, S=2, multiplicity 3(k+1), k>=1 here
    for k in range(1, nmax + 1):
        add_mode(k, 2, 3 * (k + 1), False)
    # chiral fermion: p=k+1, S=1, multiplicity 3(k+1)
    for k in range(0, nmax):
        add_mode(k + 1, 1, 3 * (k + 1), True)
    # vector boson: p=k+2, S=0, multiplicity k+1
    for k in range(0, max(0, nmax - 1)):
        add_mode(k + 2, 0, k + 1, False)
    # vector fermion: original k>=2 -> p=k-1, S=3, multiplicity k+1=p+2
    for p in range(1, nmax + 1):
        add_mode(p, 3, p + 2, True)
    return D


@njit(cache=True)
def _recurrence_one_theta(D, z0):
    nmax = D.shape[0] - 1
    Smax = D.shape[1] - 1
    Z = np.zeros((nmax + 1, Smax + 1), dtype=np.complex128)
    Z[0, :] = z0
    for n in range(1, nmax + 1):
        for k in range(1, n + 1):
            for ds in range(Smax + 1):
                d = D[k, ds]
                if d == 0:
                    continue
                for s in range(ds, Smax + 1):
                    Z[n, s] += d * Z[n - k, s - ds]
        Z[n, :] /= n
    return Z


def compute_rectangular_strip(nmax: int, Qmax: float, n_theta: int = 64):
    """Approximate Haar-projected coefficients c(n,S) for n=2J and S=3Q.

    This is a Fourier-quadrature reference implementation.  Production use
    should rerun with increasing n_theta and retain digits stable under that
    refinement.
    """
    Smax = int(round(3 * Qmax))
    out = np.zeros((nmax + 1, Smax + 1), dtype=np.float64)
    for it in range(n_theta):
        theta = 2.0 * math.pi * (it + 0.5) / n_theta
        D = build_D_theta(nmax, Smax, theta)
        z0 = zero_spin_base(Smax, theta)
        Z = _recurrence_one_theta(D, z0)
        # Relative U(2) Haar measure: (1-cos theta) dtheta/(2pi)
        out += (1.0 - math.cos(theta)) * Z.real / n_theta
    return out
