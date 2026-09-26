#!/usr/bin/env python3
"""Direct fixed-Q Appendix-A crossing test."""

import argparse
from fractions import Fraction
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from u2bps.crossings import crossing_for_sector

p=argparse.ArgumentParser()
p.add_argument("--input", type=Path, default=Path("data/reference/u2_coefficients_q_to_50.csv"))
p.add_argument("--Q", default="4/3", help="fixed microscopic charge, e.g. 4/3 or 3")
p.add_argument("--outdir", type=Path, default=Path("results/generated/appendixA_direct"))
a=p.parse_args(); a.outdir.mkdir(parents=True,exist_ok=True)
Q=float(Fraction(a.Q))
df=pd.read_csv(a.input,dtype={"coefficient":str})
# normalize exact q50 column names
if "Q" not in df:
    df["Q"]=df["y_numerator_over_3"]/3
    df["j"]=df["q_numerator_over_3"]/3
    df["twoJ"]=df["j"]-df["Q"]
res=crossing_for_sector(df,Q)
print(res)

g=df[abs(df.Q-Q)<1e-12].copy().sort_values("twoJ")
import math
from u2bps.saddle import log_dsp
g["logd"]=g.coefficient.map(lambda x:math.log(int(x)))
g["logdsp"]=(g.twoJ+Q).map(log_dsp)
plt.figure(figsize=(8,5)); plt.plot(g.twoJ,g.logd,label="exact log d"); plt.plot(g.twoJ,g.logdsp,label="single-saddle log d_sp")
if res: plt.axvline(res["twoJ_cross"],linestyle=":",label=f"crossing {res['twoJ_cross']:.4f}")
plt.xlabel("2J"); plt.ylabel("log count"); plt.legend(); plt.tight_layout();
plt.savefig(a.outdir/f"Q_{a.Q.replace('/','over')}.png",dpi=180); plt.close()
