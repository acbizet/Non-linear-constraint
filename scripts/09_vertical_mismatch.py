#!/usr/bin/env python3
"""Vertical-mismatch / shrinking-slope diagnostic.

For each analytic nonlinear-constraint point measure

 Delta_log = log d - log d_sp,
 s = d/d(2J)(log d-log d_sp),

and verify delta(2J) ~= -Delta_log/s.
"""

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from u2bps.saddle import log_dsp

p=argparse.ArgumentParser(); p.add_argument("--strip",type=Path,default=Path("data/reference/u2_thirds_Q18_J2000_clean.csv")); p.add_argument("--crossings",type=Path,default=Path("data/reference/u2_crossings_Q_thirds_to_40.csv")); p.add_argument("--outdir",type=Path,default=Path("results/generated/vertical_mismatch")); a=p.parse_args(); a.outdir.mkdir(parents=True,exist_ok=True)
raw=pd.read_csv(a.strip).rename(columns={"n":"twoJ"}); cross=pd.read_csv(a.crossings); rows=[]
for _,r in cross[cross.Q<=18+1e-12].iterrows():
    Q=float(r.Q); S=int(round(3*Q)); Jc=float(r.twoJ_ctr); Jx=float(r.twoJ_cross)
    g=raw[(raw.S==S)&(raw.coef>0)].copy(); g["logd"]=np.log(g.coef.astype(float)); center=int(round(Jc)); h=g[(g.twoJ>=center-6)&(g.twoJ<=center+6)]
    if len(h)<8: continue
    spl=CubicSpline(h.twoJ,h.logd); ld=float(spl(Jc)); sd=float(spl(Jc,1)); eps=1e-3; lsp=log_dsp(Jc+Q); ssp=(log_dsp(Jc+Q+eps)-log_dsp(Jc+Q-eps))/(2*eps); Delta=ld-lsp; s=sd-ssp
    rows.append((Q,r.Qctr,Jc,Jx,Delta,s,Jx-Jc,-Delta/s))
d=pd.DataFrame(rows,columns=["Q","Qctr","twoJ_ctr","twoJ_cross","Delta_log","slope_difference","delta_observed","delta_linear"]); d.to_csv(a.outdir/"vertical_mismatch.csv",index=False)
sub=d[d.Q>=3]; pwr=np.polyfit(np.log(sub.Qctr),np.log(abs(sub.slope_difference)),1)
plt.figure(figsize=(9,5)); plt.plot(sub.Qctr,sub.Delta_log,marker="o",markersize=2); plt.axhline(0,ls="--"); plt.xlabel("Qctr"); plt.ylabel("Delta_log"); plt.tight_layout(); plt.savefig(a.outdir/"Delta_log.png",dpi=180); plt.close()
plt.figure(figsize=(9,5)); plt.loglog(sub.Qctr,abs(sub.slope_difference),marker="o",ls="none",label="data"); x=np.linspace(sub.Qctr.min(),sub.Qctr.max(),300); C=np.exp(pwr[1]); plt.loglog(x,C*x**pwr[0],ls="--",label=fr"fit $R^{{{pwr[0]:.3f}}}$"); plt.legend(); plt.tight_layout(); plt.savefig(a.outdir/"slope_scaling.png",dpi=180); plt.close()
print(f"slope exponent = {pwr[0]:.4f}")
