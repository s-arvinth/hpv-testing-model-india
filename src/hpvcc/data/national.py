"""
Inputs of the national (India) model.

* Women by age band: NDMC India projections 2012–2070, extrapolated back to 2002 by a
  linear fit to 2012–2019.
* Rate denominator: NDMC all-age women (sum over states) for 2012–2036, extended to 2070
  at the 2031–2036 compound growth rate and back-cast to 2002 linearly.
* Recruitment Λ(t): the 18–29 NDMC series shifted ten years (Λ(2002) = 18–29 women in
  2012), extended beyond 2070 at its last five-year growth rate; the value for calendar
  year ⌊t⌋ is used (step function).
* Initial conditions: I = 20.9 % of all-age women (linear back-cast of 2012–2019 to 2002);
  C = 800 thousand prevalent cancers; C_D = 0.
"""
from __future__ import annotations

import numpy as np

from ..constants import CALIB_END, CALIB_START, SIM_END, SIM_START
from ..inputs import RecruitmentInflow, UnitInputs
from ..units import INDIA
from .sources import RawData, linear_backcast

#: Initial HPV prevalence among all-age women (%).
NATIONAL_HPV_PREVALENCE_PCT = 20.9
#: Prevalent cervical cancers in 2002 (thousands).
NATIONAL_INITIAL_CANCERS = 800.0

CALIB_YEARS = np.arange(CALIB_START, CALIB_END + 1)
_BACK_YEARS = np.arange(SIM_START, CALIB_START)


def _denominator(raw: RawData) -> dict[int, float]:
    """All-age women by year (thousands), 2002–2070."""
    total = raw.ndmc_statewise_female.groupby("Year")["Total"].sum() / 1000.0
    growth = (total[2036] / total[2031]) ** (1 / 5) - 1
    N = {y: total[y] for y in range(2012, 2037)}
    for y in range(2037, SIM_END + 1):
        N[y] = N[y - 1] * (1 + growth)
    fit = np.array([N[y] for y in CALIB_YEARS])
    from scipy.stats import linregress
    slope, intercept, *_ = linregress(CALIB_YEARS, fit)
    for y in _BACK_YEARS:
        N[y] = max(slope * y + intercept, 0.0)
    return N


def build_national_inputs(raw: RawData | None = None) -> UnitInputs:
    """Build the :class:`~hpvcc.inputs.UnitInputs` of the national model."""
    raw = raw or RawData()
    n_calib = len(CALIB_YEARS)

    # NDMC age bands, 2012-2070, back-cast to 2002.
    nd = raw.ndmc_india_bands
    i0 = int(np.where(nd["Year"].values == CALIB_START)[0][0])
    i_end = int(np.searchsorted(nd["Year"].values, SIM_END, side="right"))
    series = {k: nd[c].to_numpy(float) for k, c in (("S18", "Pop_18-29"), ("S30", "Pop_30-65"), ("S61", "Pop_66-75+"))}
    pop = {k: np.concatenate([linear_backcast(CALIB_YEARS, v[i0:i0 + n_calib], _BACK_YEARS), v[i0:i_end].copy()])
           for k, v in series.items()}

    # Recruitment: 18-29 series 2012-2070 extended ten more years at its last 5-year growth.
    s18 = series["S18"][i0:i_end].copy()
    growth = (s18[-1] / s18[-6]) ** (1 / 5) - 1
    s18_ext = np.concatenate([s18, np.array([s18[-1] * (1 + growth) ** k for k in range(1, 11)])])
    inflow_years = np.arange(SIM_START, SIM_START + len(s18_ext), dtype=float)

    N = _denominator(raw)
    denominator = np.array([N[y] for y in range(SIM_START, SIM_END + 1)])

    # Initial infections: 20.9% of all-age women, back-cast to 2002.
    infected = NATIONAL_HPV_PREVALENCE_PCT / 100.0 * np.array([N[y] for y in CALIB_YEARS])
    infected_series = np.concatenate([linear_backcast(CALIB_YEARS, infected, _BACK_YEARS), infected])
    I0 = infected_series[0]

    years = np.arange(SIM_START, SIM_END + 1, dtype=float)
    y0 = np.array([pop["S18"][0], pop["S30"][0], pop["S61"][0], 0.0, 0.0,
                   I0, NATIONAL_INITIAL_CANCERS, 0.0, 0.0, 0.0], dtype=float)

    provenance = {
        "population": f"{RawData.NDMC_INDIA_BANDS} (2012-2070), linear back-cast to {SIM_START}",
        "denominator": f"all-age women, {RawData.NDMC_STATEWISE} (2012-2036), CAGR 2031-36 to {SIM_END}",
        "recruitment": "NDMC 18-29 series shifted 10 years, step function",
        "I0": f"{NATIONAL_HPV_PREVALENCE_PCT}% of all-age women, back-cast to {SIM_START}",
        "C0": f"{NATIONAL_INITIAL_CANCERS} thousand",
    }
    return UnitInputs(INDIA, years, pop, denominator, RecruitmentInflow(inflow_years, s18_ext, "step"),
                      y0, provenance)
