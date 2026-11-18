"""
Inputs of the state/UT models.

* Women by age band and the recruitment inflow come from the state population
  projections (``popproj_female_agebands_states.xlsx``), linearly interpolated to the
  simulation years. Recruitment Λ(t) is the 18–29 population ten years ahead.
* Initial undetected infections I(2002) = literature HPV prevalence × NDMC women in the
  three age bands, the NDMC series being extrapolated back from 2012 to 2002 by a linear
  fit to 2012–2019.
* Initial cancers C(2002) = IHME GBD prevalent cases in 2002 (grouped IHME locations are
  split by fixed population shares).
* The combined UT (Dadra and Nagar Haveli and Daman and Diu) sums the initial
  conditions of its two former UTs.
"""
from __future__ import annotations

import numpy as np

from ..constants import CALIB_END, CALIB_START, SIM_END, SIM_START
from ..inputs import RecruitmentInflow, UnitInputs
from ..units import (COMBINED_UT, COMBINED_UT_PARTS, IHME_GROUP_SHARES, MODEL_TO_IHME,
                     MODEL_TO_NDMC, MODEL_TO_POPPROJ)
from .sources import RawData, linear_backcast

#: Recruitment into the 18–29 band at year t = the 18–29 population at t + this lead.
RECRUITMENT_LEAD_YEARS = 10

#: India's prevalent cervical cancers in 2012 (thousands) — fallback when IHME lacks a unit.
INDIA_PREVALENCE_2012 = 395.796

SIM_YEARS = np.arange(SIM_START, SIM_END + 1)
CALIB_YEARS = np.arange(CALIB_START, CALIB_END + 1)
_BACK_YEARS = np.arange(SIM_START, CALIB_START)
_FIT_YEARS = np.arange(CALIB_START, CALIB_END + 1)


def _ihme_prevalence(raw: RawData, unit: str, year: int) -> float | None:
    """Prevalent cancers (thousands) of ``unit`` in ``year`` from IHME GBD."""
    loc = MODEL_TO_IHME.get(unit, unit)
    value = raw.ihme_prevalence.get(loc, {}).get(int(np.clip(year, 1990, 2023)))
    if value is None:
        return None
    group = IHME_GROUP_SHARES.get(loc, {})
    if group:
        return value * (group.get(unit, 0.0) / sum(group.values())) / 1000.0
    return value / 1000.0


def _initial_conditions(raw: RawData, unit: str) -> tuple[float, float, float]:
    """(I, C, C_D) in 2002 for a single (not combined) unit, thousands."""
    hpv = raw.initial_hpv.loc[unit]
    share = float(hpv["population_share_pct"]) / 100.0
    hpv_frac = float(hpv["hpv_prevalence_pct"]) / 100.0

    bands = raw.ndmc_state_bands
    rows = bands[bands["State"] == MODEL_TO_NDMC.get(unit, unit)].sort_values("Year").reset_index(drop=True)
    if rows.empty or CALIB_START not in rows["Year"].values:
        india = raw.ndmc_india_bands_from_state_file
        i0 = int(np.where(india["Year"].values == CALIB_START)[0][0])
        cols = {k: share * india[c].to_numpy(float)[i0:] / 1000.0
                for k, c in (("S18", "Pop_18-29"), ("S30", "Pop_30-65"), ("S61", "Pop_66-75+"))}
    else:
        i0 = int(np.where(rows["Year"].values == CALIB_START)[0][0])
        cols = {k: rows[c].to_numpy(float)[i0:] / 1000.0
                for k, c in (("S18", "Pop_18-29"), ("S30", "Pop_30-65"), ("S61", "Pop_66-75+"))}
    n_calib = len(CALIB_YEARS)
    full = {k: np.concatenate([linear_backcast(_FIT_YEARS, v[:n_calib], _BACK_YEARS), v]) for k, v in cols.items()}
    infected0 = (hpv_frac * (full["S18"] + full["S30"] + full["S61"]))[0]

    cancer0 = _ihme_prevalence(raw, unit, SIM_START)
    if cancer0 is None:
        cancer0 = share * INDIA_PREVALENCE_2012
    return float(infected0), float(cancer0), 0.0


def _combined_ut(raw: RawData) -> tuple[float, float, float]:
    """Initial conditions of the combined UT: the sum of its two former UTs."""
    parts = [_initial_conditions(raw, p) for p in COMBINED_UT_PARTS]
    return tuple(a + b for a, b in zip(*parts))


def build_state_inputs(unit: str, raw: RawData | None = None) -> UnitInputs:
    """Build the :class:`~hpvcc.inputs.UnitInputs` of one state/UT model."""
    raw = raw or RawData()
    I0, C0, CD0 = _combined_ut(raw) if unit == COMBINED_UT else _initial_conditions(raw, unit)

    # Population projections: women by band, recruitment = 18–29 population ten years ahead.
    pp = raw.population_projections
    rows = pp[pp["State"] == MODEL_TO_POPPROJ.get(unit, unit)].set_index("Year").sort_index()
    if rows.empty:
        raise KeyError(f"no population projection rows for {unit!r}")
    proj_years = rows.index.to_numpy().astype(float)
    bands = {k: rows[c].to_numpy() / 1000.0 for k, c in (("S18", "18-29"), ("S30", "30-65"), ("S61", "66+"))}
    years = SIM_YEARS.astype(float)
    pop = {k: np.interp(years, proj_years, v) for k, v in bands.items()}
    recruitment = np.interp(years + RECRUITMENT_LEAD_YEARS, proj_years, bands["S18"])
    denominator = pop["S18"] + pop["S30"] + pop["S61"]
    y0 = np.array([pop["S18"][0], pop["S30"][0], pop["S61"][0], 0.0, 0.0, I0, C0, 0.0, 0.0, CD0])

    provenance = {
        "population": f"{RawData.POP_PROJECTIONS} (state '{MODEL_TO_POPPROJ.get(unit, unit)}')",
        "recruitment": "18-29 population ten years ahead, linear interpolation",
        "denominator": "women aged 18-85 (sum of the three model bands)",
        "I0": f"literature HPV prevalence ({RawData.HPV_PREVALENCE}) x NDMC women, back-cast to {SIM_START}",
        "C0": f"IHME GBD prevalent cases in {SIM_START} ({RawData.IHME})",
    }
    return UnitInputs(unit, years, pop, denominator, RecruitmentInflow(years, recruitment, "linear"),
                      y0, provenance)
