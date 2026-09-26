#!/usr/bin/env python3
"""Test the phase of the simple final-omega Gaussian prefactor.

This is deliberately labelled a diagnostic: it is not the full matrix-integral
one-loop determinant.
"""

import argparse, numpy as np
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from u2bps.saddle import saddle_data, twoJ_constraint

p=argparse.ArgumentParser(); p.add_argument("--input",type=Path,default=Path("data/reference/u2_crossings_Q_thirds_to_40.csv")); p.add_argument("--outdir",type=Path,default=Path("results/generated/gaussian")); a=p.parse_args(); a.outdir.mkdir(parents=True,exist_ok=True)
d=pd.read_csv(a.input)

def qlead(j): return -saddle_data(j)["S0"].imag/np.pi
def qgauss(j):
    s=saddle_data(j); return -s["S0"].imag/np.pi + np.angle(s["Fpp"])/(2*np.pi)
def J_for_R(R,fn):
    hi=max(20,3*R**1.5+3*R+50)
    while fn(hi)<R: hi*=2
    j=brentq(lambda x:fn(x)-R,0.1,hi)
    return j-R

rows=[]
for _,r in d.iterrows():
    Q,Jx=float(r.Q),float(r.twoJ_cross)
    best0=bestg=None
    for m in range(500):
        R=Q+2*m
        J0=J_for_R(R,qlead); Jg=J_for_R(R,qgauss)
        a0=(abs(Jx-J0),m,R,J0); ag=(abs(Jx-Jg),m,R,Jg)
        if best0 is None or a0<best0: best0=a0
        if bestg is None or ag<bestg: bestg=ag
        if J0>Jx+150 and m>5: break
    rows.append((Q,Jx,best0[2],best0[3],Jx-best0[3],bestg[2],bestg[3],Jx-bestg[3]))
r=pd.DataFrame(rows,columns=["Q","twoJ_cross","leading_Qctr","leading_twoJctr","leading_delta","gaussian_Qctr","gaussian_twoJctr","gaussian_delta"]); r.to_csv(a.outdir/"comparison.csv",index=False)
plt.figure(figsize=(9,5)); plt.plot(r.Q,abs(r.leading_delta),label="leading"); plt.plot(r.Q,abs(r.gaussian_delta),label="+ final-omega Gaussian phase"); plt.xlabel("Q"); plt.ylabel("absolute 2J residual"); plt.legend(); plt.tight_layout(); plt.savefig(a.outdir/"gaussian_comparison.png",dpi=180); plt.close()
