"""
Observed crude rates the models were fitted to (calibration window 2012–2019).

The cancer registry estimates of ICMR's National Centre for Disease Informatics and
Research / National Institute of Cancer Epidemiology and Research (NCDIR / NINE), as crude
incidence and mortality rates per 100,000 women, one row per unit and year
(``data/registry/observed_rates_2012_2019.csv``). They are used only to *show* how well
each fitted model reproduces them (Figs 1, 2, S5–S10; Table 1; Table S9). Fitting is not
part of this package; every projection runs from the fitted parameters in ``params/``.

>>> from hpvcc.observed import ObservedRates, model_fit
>>> ObservedRates.load("India").cir
>>> model_fit("India")               # model vs observed, 2012–2019
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import numpy as np
import pandas as pd

from .paths import OBSERVED_RATES_FILE
from .units import resolve_unit

#: Units whose rates are directly observed by a population-based cancer registry (PBCR);
#: the rates of all other units are PBCR-based estimates.
PBCR_OBSERVED = {"Tamil Nadu", "Mizoram", "Manipur", "Sikkim", "Tripura"}

CALIBRATION_YEARS = tuple(range(2012, 2020))


def observed_label(unit: str) -> str:
    """Legend label of the observed series of ``unit``."""
    return "NINE (PBCR)" if unit in PBCR_OBSERVED else "NINE (PBCR-estimated)"


@lru_cache(maxsize=1)
def observed_table() -> pd.DataFrame:
    """Every unit's observed CIR and CMR (columns ``Unit, Year, CIR, CMR``)."""
    return pd.read_csv(OBSERVED_RATES_FILE, float_precision="round_trip")


@dataclass(frozen=True)
class ObservedRates:
    """Observed crude incidence (``cir``) and mortality (``cmr``) rates of one unit."""

    unit: str
    years: np.ndarray
    cir: np.ndarray
    cmr: np.ndarray

    @classmethod
    def load(cls, unit: str) -> "ObservedRates":
        unit = resolve_unit(unit)
        t = observed_table()
        t = t[t.Unit == unit].sort_values("Year")
        if t.empty:
            raise KeyError(f"no observed rates for {unit!r}")
        return cls(unit, t.Year.to_numpy(float), t.CIR.to_numpy(float), t.CMR.to_numpy(float))

    @property
    def label(self) -> str:
        return observed_label(self.unit)


def model_fit(unit: str, model=None) -> pd.DataFrame:
    """Fitted (no-screening) model rates next to the observed ones, 2012–2019.

    ``model`` defaults to the unit's fitted model; pass a modified
    :class:`~hpvcc.model.HPVCervicalCancerModel` to see how a parameter change moves the fit.
    """
    from .metrics import rates
    from .model import HPVCervicalCancerModel
    unit = resolve_unit(unit)
    m = model or HPVCervicalCancerModel.for_unit(unit)
    cir, cmr = rates(m.simulate(), m.params, m.inputs.denominator)
    obs = ObservedRates.load(unit)
    yrs = list(m.years[1:].astype(int))
    idx = [yrs.index(int(y)) for y in obs.years]
    return pd.DataFrame({"Unit": unit, "Year": obs.years.astype(int),
                         "Model_CIR": cir[idx], "Observed_CIR": obs.cir,
                         "Model_CMR": cmr[idx], "Observed_CMR": obs.cmr})
