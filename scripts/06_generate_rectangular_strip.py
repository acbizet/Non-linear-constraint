#!/usr/bin/env python3
"""Reference direct (2J,3Q) high-spin recurrence with gauge Fourier quadrature.

Start with modest nmax (100--300) and increase n_theta until the Haar-projected
coefficients stabilize.  The Q<=40 production data in data/reference were made
with the same variable reorganization and higher-performance banded/chunked
runs; this readable implementation is included so the algorithm is explicit.
"""

import argparse
from pathlib import Path
import pandas as pd
from u2bps.rectangular import compute_rectangular_strip

p=argparse.ArgumentParser()
p.add_argument("--nmax",type=int,default=200)
p.add_argument("--Qmax",type=float,default=6.0)
p.add_argument("--n-theta",type=int,default=64)
p.add_argument("--output",type=Path,default=Path("results/generated/rectangular_strip.csv"))
a=p.parse_args(); a.output.parent.mkdir(parents=True,exist_ok=True)
C=compute_rectangular_strip(a.nmax,a.Qmax,a.n_theta)
rows=[]
for n in range(C.shape[0]):
    for S in range(C.shape[1]):
        if abs(C[n,S])>1e-9:
            rows.append((n,S,S/3.0,C[n,S]))
pd.DataFrame(rows,columns=["twoJ","S","Q","coefficient"]).to_csv(a.output,index=False)
print(a.output)
