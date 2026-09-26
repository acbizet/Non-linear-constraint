#!/usr/bin/env python3
"""Final convergence diagnostic: 2J_cross / 2J_ctr^closest."""

import argparse
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
p=argparse.ArgumentParser(); p.add_argument("--input",type=Path,default=Path("data/reference/u2_crossings_Q_thirds_to_40.csv")); p.add_argument("--outdir",type=Path,default=Path("results/generated/ratio")); a=p.parse_args(); a.outdir.mkdir(parents=True,exist_ok=True)
d=pd.read_csv(a.input); d["ratio"]=d.twoJ_cross/d.twoJ_ctr; d.to_csv(a.outdir/"ratio.csv",index=False)
plt.figure(figsize=(10,6)); plt.plot(d.Q,d.ratio,marker="o",markersize=3,linewidth=1.1); plt.axhline(1,ls="--"); plt.ylim(0,1.5); plt.xlabel("fixed microscopic charge Q"); plt.ylabel(r"$2J_\times/2J_{ctr}^{closest}$"); plt.title("Ratio test of the nonlinear-constraint crossing"); plt.tight_layout(); plt.savefig(a.outdir/"ratio_Q40_ylim0_1p5.png",dpi=220); plt.close()
