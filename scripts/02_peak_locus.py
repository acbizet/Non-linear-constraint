#!/usr/bin/env python3
"""Reproduce the fixed-j positive-coefficient peak analysis.

This is a useful characterization of the exact BPS distribution, but it is NOT
the direct Appendix-A nonlinear-constraint test.  Appendix A fixes Q and scans
2J; the fixed-j peak fixes j and maximizes over Q.
"""

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

p=argparse.ArgumentParser()
p.add_argument("--input", type=Path, default=Path("data/reference/u2_q100_profiles_exact.csv"))
p.add_argument("--outdir", type=Path, default=Path("results/generated/peak_locus"))
a=p.parse_args(); a.outdir.mkdir(parents=True, exist_ok=True)

df=pd.read_csv(a.input, dtype={"coefficient":str})
qcol="q_power"; qcharge="y_power"
rows=[]
for j,g in df.groupby(qcol):
    i=g["coefficient"].map(int).idxmax(); r=g.loc[i]
    rows.append((j,float(r[qcharge]),int(r["coefficient"])))
mx=pd.DataFrame(rows,columns=["j","Q_peak","max_coefficient"])
hi=mx[(mx.j>=50)&(mx.j<=100)]
xpow=hi.j.to_numpy(float)**(2/3)
A=float((xpow @ hi.Q_peak.to_numpy(float)) / (xpow @ xpow))

plt.figure(figsize=(8,5))
plt.plot(mx.j,mx.Q_peak,marker="o",markersize=2,label="exact peak")
x=np.linspace(max(1,mx.j.min()),mx.j.max(),400)
plt.plot(x,A*x**(2/3),linestyle="--",label=fr"{A:.4f} $j^{{2/3}}$")
plt.xlabel("q-power exponent j"); plt.ylabel("Q at coefficient maximum")
plt.legend(); plt.tight_layout(); plt.savefig(a.outdir/"peak_locus.png",dpi=180); plt.close()

# log peak-height fit log cmax = a j^(2/3)+b log j+c on exponent range j=50..100
X=np.column_stack([hi.j**(2/3),np.log(hi.j),np.ones(len(hi))])
y=np.log(hi.max_coefficient.map(float))
coef=np.linalg.lstsq(X,y,rcond=None)[0]
print(f"Q_peak ~ {A:.6f} j^(2/3) over q-power exponents j=50..100")
print(f"log cmax ~ {coef[0]:.6f} j^(2/3) + {coef[1]:.6f} log j + {coef[2]:.6f}")
mx.to_csv(a.outdir/"maxima.csv",index=False)
