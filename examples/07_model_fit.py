"""
07 — How well the fitted models reproduce the observed rates (2012–2019).

Compares each fitted no-screening model with the cancer registry rates it was fitted to
(ICMR-NCDIR / NINE; data/registry/). The figures are the paper's Fig 1 (India) and its
Tamil Nadu panel of Fig 2 (a state).

    python examples/07_model_fit.py
"""
from pathlib import Path

import pandas as pd

from hpvcc.observed import model_fit
from hpvcc.units import all_units
from hpvcc.viz.plots import plot_fit

OUT = Path("results/examples")

print(model_fit("India").round(2).to_string(index=False))
print(f"\nfigure -> {plot_fit('India', OUT / '07_india_fit.png')}")
print(f"figure -> {plot_fit('Tamil Nadu', OUT / '07_tamil_nadu_fit.png')}")

# Mean relative difference, model vs observed, over the calibration window, every unit
rows = []
for unit in all_units():
    m = model_fit(unit).mean(numeric_only=True)
    rows.append({"Unit": unit, "CIR diff (%)": 100 * (m.Model_CIR / m.Observed_CIR - 1),
                 "CMR diff (%)": 100 * (m.Model_CMR / m.Observed_CMR - 1)})
print("\n" + pd.DataFrame(rows).round(1).to_string(index=False))
