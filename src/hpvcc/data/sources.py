"""
Readers for the raw source files in ``data/raw`` (one method per source, each cached).

======================================================  ===========================================
file                                                    content
======================================================  ===========================================
``ihme_gbd2023_cervical_cancer_prevalence.csv``         IHME GBD 2023 prevalent cervical cancers by state
``ndmc_female_model_agebands_2012-2070.xlsx``           NDMC female population in the model age bands
                                                        (states + India)
``ndmc_india_female_agebands.xlsx``                     same, India only (national model)
``ndmc_statewise_population.xlsx``                      NDMC all-age female population by state
``popproj_female_agebands_states.xlsx``                 female population projections by state and
                                                        model age band, 1991–2100
                                                        (github.com/bisect-group/-india-female-pop-projections)
``state_initial_hpv_prevalence.csv``                    literature HPV prevalence used to initialise I
                                                        (Appendix Table S2) and each unit's share of the
                                                        national population (fallback inputs)
======================================================  ===========================================

The observed registry rates the models were fitted to are not model inputs; they live in
``data/registry/`` (see :mod:`hpvcc.observed`).
"""
from __future__ import annotations

from functools import cached_property
from pathlib import Path

import numpy as np
import pandas as pd

from ..paths import DATA_RAW


class RawData:
    """Lazy access to every raw source file."""

    IHME = "ihme_gbd2023_cervical_cancer_prevalence.csv"
    NDMC_BANDS = "ndmc_female_model_agebands_2012-2070.xlsx"
    NDMC_INDIA_BANDS = "ndmc_india_female_agebands.xlsx"
    NDMC_STATEWISE = "ndmc_statewise_population.xlsx"
    POP_PROJECTIONS = "popproj_female_agebands_states.xlsx"
    HPV_PREVALENCE = "state_initial_hpv_prevalence.csv"

    def __init__(self, directory: Path | None = None):
        self.dir = Path(directory) if directory else DATA_RAW

    def path(self, name: str) -> Path:
        return self.dir / name

    # ------------------------------------------------------------------ IHME
    @cached_property
    def ihme_prevalence(self) -> dict[str, dict[int, float]]:
        """Prevalent cervical cancers (number), ``{IHME location: {year: value}}``."""
        df = pd.read_csv(self.path(self.IHME))
        out: dict[str, dict[int, float]] = {}
        for _, row in df.iterrows():
            out.setdefault(row["location_name"], {})[int(row["year"])] = float(row["val"])
        return out

    # ------------------------------------------------------------------ populations
    @cached_property
    def ndmc_state_bands(self) -> pd.DataFrame:
        df = pd.read_excel(self.path(self.NDMC_BANDS), sheet_name="States_Female_Long")
        return df.sort_values(["State", "Year"]).reset_index(drop=True)

    @cached_property
    def ndmc_india_bands_from_state_file(self) -> pd.DataFrame:
        df = pd.read_excel(self.path(self.NDMC_BANDS), sheet_name="India_Female")
        return df.sort_values("Year").reset_index(drop=True)

    @cached_property
    def ndmc_india_bands(self) -> pd.DataFrame:
        df = pd.read_excel(self.path(self.NDMC_INDIA_BANDS), sheet_name="India_Female")
        return df.sort_values("Year").reset_index(drop=True)

    @cached_property
    def ndmc_statewise_female(self) -> pd.DataFrame:
        return pd.read_excel(self.path(self.NDMC_STATEWISE), sheet_name="Female")

    @cached_property
    def population_projections(self) -> pd.DataFrame:
        return pd.read_excel(self.path(self.POP_PROJECTIONS))

    # ------------------------------------------------------------------ literature
    @cached_property
    def initial_hpv(self) -> pd.DataFrame:
        """Per-unit HPV prevalence (%) and national population share (%)."""
        return pd.read_csv(self.path(self.HPV_PREVALENCE), float_precision="round_trip").set_index("unit")


def linear_backcast(fit_years: np.ndarray, values: np.ndarray, back_years: np.ndarray) -> np.ndarray:
    """Extrapolate a linear fit of ``values`` on ``fit_years`` back to ``back_years`` (floored at 0)."""
    from scipy.stats import linregress
    slope, intercept, *_ = linregress(fit_years, values)
    return np.maximum(slope * back_years + intercept, 0.0)
