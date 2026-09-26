#!/usr/bin/env python3
"""Analyze fixed-Q crossings against periodic nonlinear-constraint images."""

import argparse
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from u2bps.saddle import nearest_periodic_image, wrapped_mod2_residual

p=argparse.ArgumentParser()
p.add_argument("--crossings", type=Path, default=Path("data/reference/u2_crossings_Q_thirds_to_40.csv"))
p.add_argument("--outdir", type=Path, default=Path("results/generated/periodic_scan"))
a=p.parse_args(); a.outdir.mkdir(parents=True,exist_ok=True)

df=pd.read_csv(a.crossings)
if "twoJ_ctr" not in df:
    rows=[]
    for _,r in df.iterrows():
        m,R,Jc,d=nearest_periodic_image(r.Q,r.twoJ_cross)
        rows.append((m,R,Jc,d,wrapped_mod2_residual(r.Q,r.twoJ_cross)))
    x=pd.DataFrame(rows,columns=["m","Qctr","twoJ_ctr","delta","wrapped_mod2_residual"])
    df=pd.concat([df.reset_index(drop=True),x],axis=1)
else:
    df["wrapped_mod2_residual"]=[wrapped_mod2_residual(Q,J) for Q,J in zip(df.Q,df.twoJ_cross)]

df["ratio"]=df.twoJ_cross/df.twoJ_ctr
plt.figure(figsize=(9,5)); plt.plot(df.Q,df.twoJ_cross,marker="o",markersize=2,label="numerical crossing"); plt.plot(df.Q,df.twoJ_ctr,marker="s",markersize=2,label="nearest periodic image"); plt.legend(); plt.xlabel("Q"); plt.ylabel("2J"); plt.tight_layout(); plt.savefig(a.outdir/"crossings.png",dpi=180); plt.close()
plt.figure(figsize=(9,5)); plt.plot(df.Q,df.wrapped_mod2_residual,marker="o",markersize=2); plt.axhline(0,ls="--"); plt.xlabel("Q"); plt.ylabel("wrapped mod-2 residual"); plt.tight_layout(); plt.savefig(a.outdir/"mod2_residual.png",dpi=180); plt.close()
df.to_csv(a.outdir/"periodic_scan.csv",index=False)
