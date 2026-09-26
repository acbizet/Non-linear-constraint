"""Single-saddle and nonlinear-constraint formulas used throughout the project.

Conventions
-----------
After p=q=x and y=e^{2 pi i alpha}, the minimally refined trace is

    Z = Tr_BPS x^(J1+J2+Q) y^Q.

On the equal-spin slice J1=J2=J we use

    j = 2J + Q.

For U(2), Appendix-A's M=1, k=-1 saddle is

    F_{-1}(omega) = -(16/27) (omega-pi i)^3 / omega^2,
    S = F_{-1}(omega) - j omega + pi i Q.

The real part of the on-shell action is independent of microscopic Q because
pi i Q is purely imaginary.
"""

from __future__ import annotations

import math
import numpy as np
from scipy.optimize import brentq

A = 16.0 / 27.0
PI = math.pi


def saddle_roots(j: float) -> np.ndarray:
    """Roots of the cubic saddle equation in omega=i*pi*x.

    (j+16/27) x^3 - (16/9) x + 32/27 = 0.
    """
    return np.roots([j + A, 0.0, -3.0 * A, 2.0 * A])


def saddle_data(j: float):
    """Return the dominant complex saddle data at fixed j.

    Returns a dict with omega, S0, Fpp and x.  The selected root maximizes
    Re S0, which is the branch used in the numerical tests.
    """
    best = None
    c = 1j * PI
    for x in saddle_roots(j):
        omega = c * x
        F = -A * (omega - c) ** 3 / omega**2
        S0 = F - omega * j
        # d^2 F/d omega^2, used for the simple one-dimensional Gaussian test.
        Fpp = 6.0 * A * c**2 * (c - omega) / omega**4
        item = dict(x=x, omega=omega, S0=S0, Fpp=Fpp)
        if best is None or S0.real > best[0]:
            best = (S0.real, item)
    return best[1]


def log_dsp(j: float) -> float:
    """Leading single-saddle log-count Re S_{-1}(j)."""
    return float(saddle_data(j)["S0"].real)


def qctr_from_a(a: float) -> float:
    """Continuous nonlinear-constraint representative Q_ctr(a), N=2."""
    return 8.0 * a / (1.0 - a) ** 2


def j_from_a(a: float) -> float:
    """Continuous total charge j(a)=2J+Q_ctr, N=2."""
    return 8.0 * a * (1.0 + a) ** 2 / (1.0 - a) ** 3


def a_from_qctr(R: float) -> float:
    """Invert Q_ctr(a)=R on the physical 0<a<1 branch."""
    if R <= 0:
        raise ValueError("R must be positive")
    return (R + 4.0 - 2.0 * math.sqrt(2.0 * R + 4.0)) / R


def twoJ_constraint(R: float) -> float:
    """Analytic nonlinear-constraint curve 2J_ctr as a function of representative R."""
    a = a_from_qctr(R)
    return j_from_a(a) - R


def nearest_periodic_image(Q: float, twoJ_cross: float, m_max: int = 1000):
    """Nearest image R=Q+2m to a measured fixed-Q crossing.

    Returns (m, R, twoJ_ctr, delta), where delta = twoJ_cross-twoJ_ctr.
    This is an a-posteriori consistency test; it does not predict m.
    """
    best = None
    for m in range(m_max + 1):
        R = Q + 2.0 * m
        if R <= 0:
            continue
        Jc = twoJ_constraint(R)
        item = (abs(twoJ_cross - Jc), m, R, Jc, twoJ_cross - Jc)
        if best is None or item[0] < best[0]:
            best = item
        if Jc > twoJ_cross + 200.0 and m > 5:
            break
    _, m, R, Jc, delta = best
    return m, R, Jc, delta


def qctr_continuous_from_twoJ(twoJ: float) -> float:
    """Invert the continuous analytic curve 2J_ctr(R)=twoJ."""
    hi = max(4.0, 2.0 * (twoJ + 1.0) ** (2.0 / 3.0) + 10.0)
    while twoJ_constraint(hi) < twoJ:
        hi *= 2.0
    return brentq(lambda R: twoJ_constraint(R) - twoJ, 1e-12, hi)


def wrapped_mod2_residual(Q: float, twoJ_cross: float) -> float:
    """wrap(Q_ctr(twoJ_cross)-Q) to [-1,1)."""
    R = qctr_continuous_from_twoJ(twoJ_cross)
    return ((R - Q + 1.0) % 2.0) - 1.0


def gaussian_charge_shift(j: float) -> float:
    """Simple 1D Gaussian prefactor phase shift arg S''/(2*pi).

    This is *not* the full matrix-integral one-loop determinant.  It is only
    the fluctuation factor of the final omega saddle and is included as a
    diagnostic because it is easy to isolate.
    """
    return float(np.angle(saddle_data(j)["Fpp"]) / (2.0 * PI))
