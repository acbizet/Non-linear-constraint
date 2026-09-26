#!/usr/bin/env python3
"""Absolute vs relative oscillation scaling through the Q<=40 crossing table."""

import argparse, math
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from u2bps.saddle import twoJ_constraint

p=argparse.ArgumentParser(); p.add_argument("--input",type=Path,default=Path("data/reference/u2_crossings_Q_thirds_to_40.csv")); p.add_argument("--outdir",type=Path,default=Path("results/generated/oscillations")); a=p.parse_args(); a.outdir.mkdir(parents=True,exist_ok=True)
d=pd.read_csv(a.input); d["abs_delta"]=abs(d.twoJ_cross-d.twoJ_ctr); d["relpct"]=100*(d.twoJ_cross/d.twoJ_ctr-1)
# local peaks on the 1/3-spaced Q grid
y=d.abs_delta.to_numpy(); q=d.Q.to_numpy(); idx=[i for i in range(1,len(y)-1) if y[i]>=y[i-1] and y[i]>=y[i+1]]; pk=d.iloc[idx].copy()
pk.to_csv(a.outdir/"absolute_residual_local_peaks.csv",index=False)
plt.figure(figsize=(9,5)); plt.plot(d.Q,d.abs_delta,marker="o",markersize=2); plt.plot(pk.Q,pk.abs_delta,marker="s",label="local maxima"); plt.xlabel("Q"); plt.ylabel(r"$|2J_\times-2J_{ctr}|$"); plt.legend(); plt.tight_layout(); plt.savefig(a.outdir/"absolute_oscillations.png",dpi=180); plt.close()
# spacing-normalized residual
def spacing(R): return 0.5*((twoJ_constraint(R+2)-twoJ_constraint(R))+(twoJ_constraint(R)-twoJ_constraint(R-2)))
d["fraction_half_spacing"]=2*(d.twoJ_cross-d.twoJ_ctr)/d.Qctr.map(spacing)
plt.figure(figsize=(9,5)); plt.plot(d.Q,d.fraction_half_spacing,marker="o",markersize=2); plt.axhline(0,ls="--"); plt.axhline(1,ls=":"); plt.axhline(-1,ls=":"); plt.xlabel("Q"); plt.ylabel("2 Delta(2J) / local image spacing"); plt.tight_layout(); plt.savefig(a.outdir/"spacing_normalized.png",dpi=180); plt.close()
print(pk[["Q","Qctr","abs_delta","relpct"]].to_string(index=False))
