#!/usr/bin/env python3
"""Raw y=-1 index control experiment.

At fixed j, the exact index is sum_Q (-1)^Q d(j,Q) in integer-Q sectors.
This experiment demonstrates the enormous cancellation and explains why the raw
index magnitude is not the single-saddle envelope used in Appendix A.
"""

import argparse, math
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

p=argparse.ArgumentParser()
p.add_argument("--input", type=Path, default=Path("data/reference/u2_q100_profiles_exact.csv"))
p.add_argument("--outdir", type=Path, default=Path("results/generated/raw_index"))
a=p.parse_args(); a.outdir.mkdir(parents=True,exist_ok=True)

df=pd.read_csv(a.input,dtype={"coefficient":str})
rows=[]
for j,g in df.groupby("q_power"):
    terms=[((-1)**int(Q))*int(c) for Q,c in zip(g.y_power,g.coefficient)]
    idx=sum(terms); total=sum(int(c) for c in g.coefficient)
    rows.append((j,idx,total,abs(idx)/total if total else 0.0))
r=pd.DataFrame(rows,columns=["j","index","positive_sum","cancellation_ratio"])
r.to_csv(a.outdir/"raw_index.csv",index=False)
plt.figure(figsize=(8,5)); plt.semilogy(r.j,r.cancellation_ratio)
plt.xlabel("j"); plt.ylabel(r"$|index|/\sum_Q d(j,Q)$"); plt.tight_layout()
plt.savefig(a.outdir/"index_cancellation.png",dpi=180); plt.close()
print(r.tail(10).to_string(index=False))
