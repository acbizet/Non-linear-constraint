"""Fixed-charge crossing extraction and periodic-image diagnostics."""

from __future__ import annotations

import math
import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline

from .saddle import (
    log_dsp,
    nearest_periodic_image,
    qctr_continuous_from_twoJ,
    wrapped_mod2_residual,
)


def add_fixed_Q_coordinates(df: pd.DataFrame) -> pd.DataFrame:
    """Normalize a coefficient table to columns Q, twoJ, j, coefficient."""
    x = df.copy()
    if "Q" not in x:
        x["Q"] = x["y_numerator_over_3"] / 3.0
    if "j" not in x:
        if "q_power" in x:
            x["j"] = x["q_power"]
        else:
            x["j"] = x["q_numerator_over_3"] / 3.0
    if "twoJ" not in x:
        x["twoJ"] = x["j"] - x["Q"]
    return x


def crossing_for_sector(g: pd.DataFrame, Q: float, coefficient_col: str = "coefficient"):
    """Largest-spin sign-changing crossing of log d_sp - log d at fixed Q."""
    g = add_fixed_Q_coordinates(g)
    g = g[np.isclose(g["Q"], Q)].sort_values("twoJ").copy()
    g = g[g[coefficient_col].map(lambda v: int(v) > 0)].copy()
    if g.empty:
        return None
    g["logd"] = g[coefficient_col].map(lambda v: math.log(int(v)))
    g["logdsp"] = (g["twoJ"] + Q).map(log_dsp)
    g["diff"] = g["logdsp"] - g["logd"]

    x = g["twoJ"].to_numpy(float)
    f = g["diff"].to_numpy(float)
    roots = []
    for i in range(len(x) - 1):
        if f[i] == 0:
            roots.append(x[i])
        elif f[i] * f[i + 1] < 0:
            roots.append(x[i] - f[i] * (x[i + 1] - x[i]) / (f[i + 1] - f[i]))
    if not roots:
        return None

    Jx = max(roots)
    m, R, Jc, delta = nearest_periodic_image(Q, Jx)
    Rcont = qctr_continuous_from_twoJ(Jx)
    return {
        "Q": Q,
        "twoJ_cross": Jx,
        "j_cross": Jx + Q,
        "m": m,
        "Qctr": R,
        "twoJ_ctr": Jc,
        "delta": delta,
        "abs_delta": abs(delta),
        "relative_error_percent": 100.0 * delta / Jc,
        "Qctr_continuous": Rcont,
        "wrapped_mod2_residual": wrapped_mod2_residual(Q, Jx),
    }


def vertical_mismatch_and_slope(g: pd.DataFrame, Q: float, twoJ_ctr: float, coefficient_col="coefficient"):
    """Measure Delta_log and d/d(2J)(log d-log d_sp) at an analytic locus.

    A local cubic spline is used only to interpolate the discrete exact coefficient
    data.  This is the diagnostic used to distinguish an O(1)/logarithmic vertical
    correction from the apparently power-like horizontal crossing displacement.
    """
    g = add_fixed_Q_coordinates(g)
    g = g[np.isclose(g["Q"], Q)].sort_values("twoJ").copy()
    g = g[g[coefficient_col].map(lambda v: float(v) > 0)].copy()
    g["logd"] = g[coefficient_col].map(lambda v: math.log(float(v)))

    center = int(round(twoJ_ctr))
    h = g[(g["twoJ"] >= center - 6) & (g["twoJ"] <= center + 6)].copy()
    if len(h) < 8:
        raise ValueError("not enough points around analytic locus")
    spl = CubicSpline(h["twoJ"].to_numpy(float), h["logd"].to_numpy(float))
    logd = float(spl(twoJ_ctr))
    slope_d = float(spl(twoJ_ctr, 1))

    eps = 1e-3
    logsp = log_dsp(twoJ_ctr + Q)
    slope_sp = (log_dsp(twoJ_ctr + Q + eps) - log_dsp(twoJ_ctr + Q - eps)) / (2 * eps)
    Delta = logd - logsp
    slope = slope_d - slope_sp
    return {
        "Delta_log": Delta,
        "slope_logd": slope_d,
        "slope_logdsp": slope_sp,
        "slope_difference": slope,
        "delta_twoJ_linear": -Delta / slope,
    }
