"""Exact sparse Taylor expansion of the positive U(2) BPS counting function.

This is the transparent Python-integer implementation used to validate the
coefficient formulas.  It is practical through the coefficient of q^50
(q-power exponent j~50) on a typical workstation; reaching q^100 (j~100) is
possible but substantially heavier.

Integerized variables:

    x = r^3,  y = s^3,

so a coefficient is labelled by R=3j and S=3Q.  The U(2) adjoint weights in
the relative gauge fugacity z are {0,0,+1,-1}.

The oscillator families are

 scalar bosons:   r^(2+3k) s^2, multiplicity 3(k+1), k>=0
 chiral fermions: r^(4+3k) s,   multiplicity 3(k+1), k>=0
 vector bosons:   r^(6+3k),     multiplicity k+1,    k>=0
 vector fermions: r^(3k) s^3,   multiplicity k+1,    k>=1

Bosons contribute (1-X)^(-m), fermions (1+X)^m.  Rather than multiply the
product directly, we use the logarithmic-derivative recurrence

    R Z_R = sum_{k=1}^R D_k Z_{R-k},
    D = r d/dr log Z.

After the relative-gauge expansion, Haar projection is simply

    c_{R,S} = [z^0] Z_{R,S} - [z^1] Z_{R,S}.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
import pandas as pd

WEIGHTS = (0, 0, 1, -1)


@dataclass(frozen=True)
class Oscillator:
    r_power: int
    s_power: int
    multiplicity: int
    fermion: bool


def oscillator_families(Rmax: int):
    k = 0
    while 2 + 3 * k <= Rmax:
        yield Oscillator(2 + 3 * k, 2, 3 * (k + 1), False)
        k += 1
    k = 0
    while 4 + 3 * k <= Rmax:
        yield Oscillator(4 + 3 * k, 1, 3 * (k + 1), True)
        k += 1
    k = 0
    while 6 + 3 * k <= Rmax:
        yield Oscillator(6 + 3 * k, 0, k + 1, False)
        k += 1
    k = 1
    while 3 * k <= Rmax:
        yield Oscillator(3 * k, 3, k + 1, True)
        k += 1


def build_log_derivative(Rmax: int, Smax: int):
    """D[R][(S,z_power)] for r d/dr log Z."""
    D = [defaultdict(int) for _ in range(Rmax + 1)]
    for osc in oscillator_families(Rmax):
        nhmax = Rmax // osc.r_power
        if osc.s_power:
            nhmax = min(nhmax, Smax // osc.s_power)
        for n in range(1, nhmax + 1):
            sign = (-1) ** (n + 1) if osc.fermion else 1
            coeff = osc.multiplicity * osc.r_power * sign
            R = osc.r_power * n
            S = osc.s_power * n
            for w in WEIGHTS:
                D[R][(S, w * n)] += coeff
    return D


def compute_exact_coefficients(qmax: int, Qmax: float | None = None, progress: bool = False):
    """Compute exact positive U(2) coefficients with Python arbitrary integers.

    Parameters
    ----------
    qmax:
        Maximum exponent j of q/x in the Taylor expansion (retain terms through q^j).  Internally Rmax=3*qmax.
    Qmax:
        Optional physical-charge cutoff.  Defaults to Qmax=qmax.
    progress:
        Print the number of sparse Laurent terms at each R.
    """
    Rmax = 3 * int(qmax)
    Smax = Rmax if Qmax is None else int(round(3 * Qmax))
    D = build_log_derivative(Rmax, Smax)
    Z = [dict() for _ in range(Rmax + 1)]
    Z[0] = {(0, 0): 1}

    for R in range(1, Rmax + 1):
        acc = defaultdict(int)
        for k in range(1, R + 1):
            if not D[k] or not Z[R - k]:
                continue
            for (Sd, Zd), cd in D[k].items():
                for (Sz, Zz), cz in Z[R - k].items():
                    S = Sd + Sz
                    if S <= Smax:
                        acc[(S, Zd + Zz)] += cd * cz
        zr = {}
        for key, val in acc.items():
            q, rem = divmod(val, R)
            if rem:
                raise ArithmeticError(f"recurrence failed exact division at R={R}, key={key}")
            if q:
                zr[key] = q
        Z[R] = zr
        if progress:
            print(f"R={R:4d}/{Rmax}: {len(zr):8d} Laurent terms")

    rows = []
    for R, zr in enumerate(Z):
        for S in sorted({s for s, _ in zr.keys()}):
            c = zr.get((S, 0), 0) - zr.get((S, 1), 0)
            if c:
                rows.append(
                    {
                        "q_numerator_over_3": R,
                        "y_numerator_over_3": S,
                        "q_power": R / 3.0,
                        "Q": S / 3.0,
                        "coefficient": c,
                    }
                )
    return pd.DataFrame(rows)


def fixed_q_maxima(coeffs: pd.DataFrame) -> pd.DataFrame:
    """Peak coefficient and charge at every available q power."""
    rows = []
    for q, g in coeffs.groupby("q_power"):
        i = g["coefficient"].map(int).idxmax()
        r = g.loc[i]
        rows.append(
            {
                "q_power": q,
                "Q_peak": r["Q"],
                "max_coefficient": int(r["coefficient"]),
                "sum_coefficients": sum(map(int, g["coefficient"])),
            }
        )
    return pd.DataFrame(rows).sort_values("q_power")
