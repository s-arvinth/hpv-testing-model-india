"""
06 — All 36 state/UT models side by side: CIR in 2050 without screening and at 70 %
HPV-DNA coverage (deterministic mid-range runs, a few seconds in total).

The paper's numbers come from 1000-simulation ensembles (`hpvcc paper ...`); the
deterministic runs here are close but not identical.

    python examples/06_compare_states.py
"""
import pandas as pd

from hpvcc import project
from hpvcc.units import REGION_OF, STATES

rows = []
for unit in STATES:
    s = project(unit, coverage=[0, 0.7], deterministic=True).summary(years=(2050,))
    rows.append({"Unit": unit, "Region": REGION_OF[unit],
                 "CIR 2050, no screening": s.CIR_2050.iloc[0], "CIR 2050, 70 %": s.CIR_2050.iloc[1]})
t = pd.DataFrame(rows).sort_values("CIR 2050, 70 %", ascending=False)
t["below 4 at 70 %"] = t["CIR 2050, 70 %"] < 4
print(t.round(2).to_string(index=False))
print(f"\n{t['below 4 at 70 %'].sum()} of {len(t)} units below 4 per 100,000 in 2050 at 70 % coverage")
