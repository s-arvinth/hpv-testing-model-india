r"""
Epidemiological outputs computed from trajectories.

Crude incidence rate (CIR) — *inflow method*
    New cancers in year t are the change in the cancer stock plus everything that
    left it during the year::

        cases(t) = max( C(t) − C(t−1) + (r_C + μ)·C(t−1) + μ_C·C(t−1)^n_PL , 0 )
        CIR(t)   = cases(t) / N(t) × 10⁵

    (Using ΔC alone would ignore cancers that arose and left within the year.)

Crude mortality rate (CMR)
    C_D is a cumulative counter, so deaths(t) = max(C_D(t) − C_D(t−1), 0) and
    CMR(t) = deaths(t) / N(t) × 10⁵.

Ensemble intervals
    For a stochastic ensemble the paper's figures use the *trajectory method*: take the
    mean and the 2.5th / 97.5th percentiles of the C and C_D trajectories across
    simulations, convert each to a rate, and report ``mean ± |hi − lo| / 2``.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .constants import PER_100K
from .parameters import ModelParameters


def incidence_counts(C: np.ndarray, params: ModelParameters) -> np.ndarray:
    """Annual new cancers (thousands) along the last axis of ``C`` (length n → n−1)."""
    C = np.asarray(C, dtype=float)
    outflow = (params.r_C + params.mu) * C[..., :-1] + params.mu_C * np.maximum(C[..., :-1], 0.0) ** params.n_PL
    return np.maximum(np.diff(C, axis=-1) + outflow, 0.0)


def death_counts(C_D: np.ndarray) -> np.ndarray:
    """Annual cancer deaths (thousands) from the cumulative counter ``C_D``."""
    return np.maximum(np.diff(np.asarray(C_D, dtype=float), axis=-1), 0.0)


def rates(traj: np.ndarray, params: ModelParameters, denominator: np.ndarray):
    """CIR and CMR per 100,000 for years 2..n of a trajectory of shape ``(n, 10)``."""
    N = np.maximum(np.asarray(denominator, dtype=float)[1:len(traj)], 1e-12)
    cir = incidence_counts(traj[:, 6], params) / N * PER_100K
    cmr = death_counts(traj[:, 9]) / N * PER_100K
    return cir, cmr


# --------------------------------------------------------------------------------------
# Ensemble summaries
# --------------------------------------------------------------------------------------

@dataclass
class TrajectoryStats:
    """Mean and 2.5/97.5 percentile trajectories of one compartment over simulations."""

    mean: np.ndarray
    low: np.ndarray
    high: np.ndarray

    @classmethod
    def of(cls, sims: np.ndarray, axis: int = 0) -> "TrajectoryStats":
        return cls(np.mean(sims, axis=axis),
                   np.percentile(sims, 2.5, axis=axis),
                   np.percentile(sims, 97.5, axis=axis))


def _inflow_scalar(c_prev: float, c_curr: float, params: ModelParameters) -> float:
    out = (params.r_C + params.mu) * c_prev + params.mu_C * max(c_prev, 0.0) ** params.n_PL
    return max((c_curr - c_prev) + out, 0.0)


@dataclass
class RateSnapshot:
    """Rates in one year with the figure convention for the 95% interval."""

    cir: float
    cir_halfwidth: float
    cmr: float
    cmr_halfwidth: float

    def as_tuple(self):
        return self.cir, self.cir_halfwidth, self.cmr, self.cmr_halfwidth


def snapshot(C: TrajectoryStats, C_D: TrajectoryStats, years: np.ndarray, year: int,
             params: ModelParameters, denominator: np.ndarray) -> RateSnapshot:
    """CIR/CMR in ``year`` from trajectory statistics (mean ± half-width of the percentile rates)."""
    iy = int(np.where(np.asarray(years) == year)[0][0])
    N = max(float(denominator[iy]), 1e-12)
    dCm = _inflow_scalar(C.mean[iy - 1], C.mean[iy], params)
    dCl = _inflow_scalar(C.low[iy - 1], C.low[iy], params)
    dCh = _inflow_scalar(C.high[iy - 1], C.high[iy], params)
    dDm = max(C_D.mean[iy] - C_D.mean[iy - 1], 0.0)
    dDl = max(C_D.low[iy] - C_D.low[iy - 1], 0.0)
    dDh = max(C_D.high[iy] - C_D.high[iy - 1], 0.0)
    return RateSnapshot(dCm / N * PER_100K, abs(dCh - dCl) / 2 / N * PER_100K,
                        dDm / N * PER_100K, abs(dDh - dDl) / 2 / N * PER_100K)
