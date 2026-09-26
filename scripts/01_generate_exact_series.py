#!/usr/bin/env python3
"""Generate exact positive U(2) coefficients with arbitrary-size integers."""

import argparse
from pathlib import Path
from u2bps.exact_series import compute_exact_coefficients, fixed_q_maxima

p = argparse.ArgumentParser()
p.add_argument("--qmax", type=int, default=50, help="maximum q-power exponent j; retains coefficients through q^j")
p.add_argument("--Qmax", type=float, default=None)
p.add_argument("--outdir", type=Path, default=Path("results/generated/exact"))
p.add_argument("--progress", action="store_true")
a = p.parse_args()
a.outdir.mkdir(parents=True, exist_ok=True)

df = compute_exact_coefficients(a.qmax, a.Qmax, progress=a.progress)
df.to_csv(a.outdir / f"u2_coefficients_q_to_{a.qmax}.csv", index=False)
mx = fixed_q_maxima(df)
mx.to_csv(a.outdir / f"u2_maxima_q_to_{a.qmax}.csv", index=False)
print(mx.tail(10).to_string(index=False))
