"""
Model fit figures: fitted no-screening model vs the observed registry rates, 2012–2019,
drawn by :func:`hpvcc.viz.figures.fit_bars`.

* Fig 1 — India (CIR and CMR side by side)
* Fig 2 and Figs S5–S10 — one state/UT per panel
"""
from __future__ import annotations

from pathlib import Path

from ...observed import model_fit, observed_label
from ...viz.figures import fit_bars


def fit_figure(unit: str, out: Path) -> Path:
    f = model_fit(unit)
    return fit_bars(unit, f.Year, f.Model_CIR.to_numpy(), f.Observed_CIR.to_numpy(), f.Model_CMR.to_numpy(),
                    f.Observed_CMR.to_numpy(), out, observed_label=observed_label(unit))
