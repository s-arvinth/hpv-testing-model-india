"""
Master results tables, built straight from the cached ensembles.

``MasterTables(...).build()`` returns the data frames written to
``results/paper/tables/csv/`` (one CSV each):

=================================  ==========================================================
``Annual_Long``                    unit × coverage × year: CIR/CMR (mean, 95 % CI), cases,
                                   deaths, reductions, averted (annual + cumulative from 2026)
``Coverage_Required``              minimum coverage with mean CIR < 4, per unit and year
``Elimination_Year``               first / sustained year with mean CIR < 4, per coverage
``Calibration_2012_2019``          fitted model vs observed registry rates, 2012–2019
``Cumulative_Summary``             cumulative cases/deaths averted 2026 → horizon
``National_Totals_from_States``    sums over the 36 state models
``Parameters``                     final parameters of every unit
``Deterministic_Baseline``         no-screening trajectories of all compartments
``India_TestType``                 national test-type comparison at 75 % coverage (Table S8)
``TestType_by_Unit``               unit × screening test × year at 70 % coverage: CIR/CMR
                                   (mean, 95 % CI) and reductions from no screening
``TestType_Elimination``           unit × screening test at 70 % coverage: first year with
                                   mean CIR < 4 and the rates in the snapshot years
=================================  ==========================================================

Interval conventions (identical to the figures): annual rates use the trajectory method
(mean ± half-width of the 2.5/97.5-percentile trajectories); cumulative quantities are
computed per simulation and summarised by mean and percentiles.
"""
from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from ..constants import SCREENING_START, SCREENING_TESTS, WHO_THRESHOLD
from .experiments import Experiment
from ..model import HPVCervicalCancerModel
from ..metrics import incidence_counts, death_counts, rates
from ..parameters import ModelParameters, load_parameters
from ..paths import TABLES_DIR, ensure
from ..observed import ObservedRates
from ..units import INDIA, REGION_OF, STATES

NATIONAL_REGION = "National (dedicated India model)"
CUMULATIVE_HORIZONS = list(range(2030, 2071, 5))
SNAPSHOT_YEARS = list(range(2030, 2071, 5))


def region(unit: str) -> str:
    return NATIONAL_REGION if unit == INDIA else REGION_OF[unit]


def _pct(a, q, axis):
    return np.percentile(a, q, axis=axis)


def unit_frames(unit: str, C: np.ndarray, CD: np.ndarray, years: np.ndarray, N: np.ndarray,
                coverages: np.ndarray, params: ModelParameters):
    """Annual, coverage-required and elimination tables of one unit from its ensemble."""
    reg = region(unit)
    inflow = lambda X: incidence_counts(X, params)            # noqa: E731
    yr = years[1:].astype(int)
    k = 1e5 / N[1:]
    covs = np.round(np.asarray(coverages) * 100).astype(int)
    Cm, Cl, Ch = C.mean(1), _pct(C, 2.5, 1), _pct(C, 97.5, 1)
    Dm, Dl, Dh = CD.mean(1), _pct(CD, 2.5, 1), _pct(CD, 97.5, 1)
    im, il, ih = inflow(Cm), inflow(Cl), inflow(Ch)
    dm, dl, dh = death_counts(Dm), death_counts(Dl), death_counts(Dh)
    hwi, hwd = np.abs(ih - il) / 2, np.abs(dh - dl) / 2
    i_s, d_s = inflow(C), death_counts(CD)
    scr = (yr >= SCREENING_START).astype(float)
    cum_c = np.cumsum(i_s * scr, axis=-1) * 1000
    cum_d = np.cumsum(d_s * scr, axis=-1) * 1000
    base_cc, base_cd = cum_c[0].mean(0), cum_d[0].mean(0)
    av_c = base_cc[None, None, :] - cum_c
    av_d = base_cd[None, None, :] - cum_d
    F = lambda a: a.reshape(-1)                                # noqa: E731
    n = len(covs) * len(yr)
    rows = {
        "Unit": np.repeat(unit, n), "Region": np.repeat(reg, n),
        "Coverage_pct": np.repeat(covs, len(yr)), "Year": np.tile(yr, len(covs)),
        "Screening_active": F(np.tile(yr >= SCREENING_START, (len(covs), 1))),
        "Women_N": F(np.tile(N[1:] * 1000, (len(covs), 1))),
        "CIR_mean": F(im * k), "CIR_lo95": F(np.maximum(im - hwi, 0) * k), "CIR_hi95": F((im + hwi) * k),
        "CIR_hw95": F(hwi * k),
        "CIR_simpct_lo95": F(_pct(i_s * k, 2.5, 1)), "CIR_simpct_hi95": F(_pct(i_s * k, 97.5, 1)),
        "CMR_mean": F(dm * k), "CMR_lo95": F(np.maximum(dm - hwd, 0) * k), "CMR_hi95": F((dm + hwd) * k),
        "CMR_hw95": F(hwd * k),
        "CMR_simpct_lo95": F(_pct(d_s * k, 2.5, 1)), "CMR_simpct_hi95": F(_pct(d_s * k, 97.5, 1)),
        "Cases_mean": F(im * 1000), "Cases_lo95": F(np.maximum(im - hwi, 0) * 1000), "Cases_hi95": F((im + hwi) * 1000),
        "Deaths_mean": F(dm * 1000), "Deaths_lo95": F(np.maximum(dm - hwd, 0) * 1000), "Deaths_hi95": F((dm + hwd) * 1000),
        "Prevalent_cancer_mean": F(Cm[:, 1:] * 1000),
        "CIR_reduction_pct_vs_baseline": F(100 * (1 - im / np.where(im[0] > 0, im[0], np.nan))),
        "CMR_reduction_pct_vs_baseline": F(100 * (1 - dm / np.where(dm[0] > 0, dm[0], np.nan))),
        "Cases_averted_annual": F((im[0][None] - im) * 1000),
        "Deaths_averted_annual": F((dm[0][None] - dm) * 1000),
        "CumCases_from2026_mean": F(cum_c.mean(1)), "CumCases_from2026_lo95": F(_pct(cum_c, 2.5, 1)),
        "CumCases_from2026_hi95": F(_pct(cum_c, 97.5, 1)),
        "CumDeaths_from2026_mean": F(cum_d.mean(1)), "CumDeaths_from2026_lo95": F(_pct(cum_d, 2.5, 1)),
        "CumDeaths_from2026_hi95": F(_pct(cum_d, 97.5, 1)),
        "CumCases_averted_mean": F(av_c.mean(1)), "CumCases_averted_lo95": F(_pct(av_c, 2.5, 1)),
        "CumCases_averted_hi95": F(_pct(av_c, 97.5, 1)),
        "CumDeaths_averted_mean": F(av_d.mean(1)), "CumDeaths_averted_lo95": F(_pct(av_d, 2.5, 1)),
        "CumDeaths_averted_hi95": F(_pct(av_d, 97.5, 1)),
        "CumCases_averted_pct": F(100 * av_c.mean(1) / np.where(base_cc > 0, base_cc, np.nan)[None]),
        "CumDeaths_averted_pct": F(100 * av_d.mean(1) / np.where(base_cd > 0, base_cd, np.nan)[None]),
    }
    long = pd.DataFrame(rows)
    long["Below_WHO_4"] = long["CIR_mean"] < WHO_THRESHOLD

    def first_cov(v):
        w = np.where(v < WHO_THRESHOLD)[0]
        return int(covs[w[0]]) if len(w) else None

    cr = []
    for t, Y in enumerate(yr):
        if Y < 2030:
            continue
        m, hw = im[:, t] * k[t], np.abs(ih[:, t] - il[:, t]) / 2 * k[t]
        req = first_cov(m)
        if req is None:
            cr.append((unit, reg, int(Y), 100, False, None, None, float(m[-1])))
        else:
            lo, hi = first_cov(m - hw), first_cov(m + hw)
            cr.append((unit, reg, int(Y), req, True, lo if lo is not None else req,
                       hi if hi is not None else 100, float(m[list(covs).index(req)])))
    covreq = pd.DataFrame(cr, columns=["Unit", "Region", "Year", "Coverage_required_pct", "Achievable_by_100pct",
                                       "Coverage_required_lo95", "Coverage_required_hi95", "CIR_at_required_coverage"])

    post = yr >= SCREENING_START

    def first_year(v):
        w = np.where((v < WHO_THRESHOLD) & post)[0]
        return int(yr[w[0]]) if len(w) else None

    def sustained(v):
        below = (v < WHO_THRESHOLD) | ~post
        if not below[-1]:
            return None
        j = len(v) - 1
        while j > 0 and below[j - 1] and post[j - 1]:
            j -= 1
        return int(yr[j]) if post[j] else SCREENING_START

    el = []
    for ai, c in enumerate(covs):
        m = im[ai] * k
        hw = np.abs(ih[ai] - il[ai]) / 2 * k
        el.append((unit, reg, int(c), first_year(m), sustained(m), first_year(m - hw), first_year(m + hw),
                   *[float(m[list(yr).index(Y)]) for Y in (2030, 2050, 2070)]))
    elim = pd.DataFrame(el, columns=["Unit", "Region", "Coverage_pct", "First_year_CIR_below4",
                                     "Sustained_below4_from", "First_year_optimistic(lo95)",
                                     "First_year_pessimistic(hi95)", "CIR_2030", "CIR_2050", "CIR_2070"])
    return long, covreq, elim


def calibration_rows(unit: str, model: HPVCervicalCancerModel, traj: np.ndarray) -> list[tuple]:
    """Fitted model vs observed registry rates for 2012–2019 from the deterministic trajectory."""
    p, inp = model.params, model.inputs
    inc = incidence_counts(traj[:, 6], p)
    dth = death_counts(traj[:, 9])
    yr = list(model.years[1:].astype(int))
    N = inp.denominator
    obs = ObservedRates.load(unit)
    out = []
    for i, Y in enumerate(obs.years.astype(int)):
        t = yr.index(Y)
        out.append((unit, int(Y), inc[t] * 1e3, dth[t] * 1e3,
                    inc[t] / N[t + 1] * 1e5, obs.cir[i], dth[t] / N[t + 1] * 1e5, obs.cmr[i]))
    return out


def parameter_row(unit: str, p: ModelParameters, n_sims: int, seed) -> dict:
    row = {"Unit": unit, "Region": region(unit), "p": p.p, "K1": p.K1, "K2": p.K2,
           "beta_const": p.beta, "beta_Reinfec": p.beta_reinfection, "gamma": p.gamma,
           "r_I": p.r_I, "r_C": p.r_C, "mu_C": p.mu_C, "n_PL": p.n_PL, "omega": p.omega,
           "epsilon": p.epsilon, "mu": p.mu, "H_pos": p.H_pos,
           "n_sims": n_sims, "seed": seed}
    return row


class MasterTables:
    """Build every master table from the cached experiments."""

    def __init__(self, states_experiment: str = "states_coverage_sweep",
                 national_experiment: str = "india_coverage_sweep",
                 test_type_experiment: str = "india_test_type",
                 states_test_type_experiment: str = "states_test_type",
                 national_test_type_experiment: str = "india_test_type_70"):
        self.states_exp = Experiment.get(states_experiment)
        self.national_exp = Experiment.get(national_experiment)
        self.test_exp = Experiment.get(test_type_experiment)
        self.states_test_exp = Experiment.get(states_test_type_experiment)
        self.national_test_exp = Experiment.get(national_test_type_experiment)

    def build(self, verbose: bool = False) -> dict[str, pd.DataFrame]:
        L, CR, EL, CAL, PAR, DET = [], [], [], [], [], []
        for unit in [*STATES, INDIA]:
            exp = self.national_exp if unit == INDIA else self.states_exp
            res = exp.result(unit)
            model = HPVCervicalCancerModel.for_unit(unit)
            cov = np.array([s.coverage for s in res.scenarios])
            l, c, e = unit_frames(unit, res.C, res.C_D, res.years, res.denominator, cov, model.params)
            L.append(l); CR.append(c); EL.append(e)
            traj = model.simulate()
            CAL += calibration_rows(unit, model, traj)
            PAR.append(parameter_row(unit, model.params, res.n_sims, res.meta.get("draws", {}).get("seed")))
            for j, Y in enumerate(model.years):
                DET.append((unit, int(Y), *traj[j], model.inputs.denominator[j]))
            if verbose:
                print(f"  tables: {unit}")
        long = pd.concat(L, ignore_index=True)
        long["CIR_rank_desc"] = long.groupby(["Year", "Coverage_pct"])["CIR_mean"].rank(ascending=False, method="min").astype(int)
        long["CMR_rank_desc"] = long.groupby(["Year", "Coverage_pct"])["CMR_mean"].rank(ascending=False, method="min").astype(int)
        with warnings.catch_warnings():   # all-None year columns for units never below 4
            warnings.simplefilter("ignore", FutureWarning)
            covreq = pd.concat(CR, ignore_index=True)
            elim = pd.concat(EL, ignore_index=True)
        calib = pd.DataFrame(CAL, columns=["Unit", "Year", "Model_cases", "Model_deaths",
                                           "Model_CIR", "Observed_CIR", "Model_CMR", "Observed_CMR"])
        det = pd.DataFrame(DET, columns=["Unit", "Year", "S18", "S30", "S61", "S_R", "H_T", "I", "C", "R_H",
                                         "C_S", "C_D_cum", "N_total"])
        cs = long[long.Year.isin(CUMULATIVE_HORIZONS)][[
            "Unit", "Region", "Coverage_pct", "Year",
            "CumCases_from2026_mean", "CumCases_from2026_lo95", "CumCases_from2026_hi95",
            "CumCases_averted_mean", "CumCases_averted_lo95", "CumCases_averted_hi95", "CumCases_averted_pct",
            "CumDeaths_from2026_mean", "CumDeaths_from2026_lo95", "CumDeaths_from2026_hi95",
            "CumDeaths_averted_mean", "CumDeaths_averted_lo95", "CumDeaths_averted_hi95", "CumDeaths_averted_pct",
        ]].rename(columns={"Year": "Horizon_year(2026..Y)"})
        st = long[long.Unit != INDIA]
        nat = st.groupby(["Year", "Coverage_pct"], as_index=False)[[
            "Women_N", "Cases_mean", "Deaths_mean", "Cases_averted_annual", "Deaths_averted_annual",
            "CumCases_from2026_mean", "CumDeaths_from2026_mean", "CumCases_averted_mean",
            "CumDeaths_averted_mean"]].sum()
        nat["CIR_popweighted"] = nat.Cases_mean / nat.Women_N * 1e5
        nat["CMR_popweighted"] = nat.Deaths_mean / nat.Women_N * 1e5
        return {"Annual_Long": long, "Coverage_Required": covreq, "Elimination_Year": elim,
                "Cumulative_Summary": cs, "National_Totals_from_States": nat,
                "Calibration_2012_2019": calib, "Parameters": pd.DataFrame(PAR),
                "Deterministic_Baseline": det, "India_TestType": self.test_type_table(),
                **self.test_type_by_unit()}

    def test_type_by_unit(self) -> dict[str, pd.DataFrame]:
        """The four screening tests at one coverage (70 %) in every unit.

        ``TestType_by_Unit`` has one row per unit, test and year (2026–2070): CIR and CMR
        (mean and 95 % interval, the convention of every other table), the no-screening
        rates of the same year and the reductions. ``TestType_Elimination`` has one row per
        unit and test: the first year with mean CIR < 4 per 100,000 and the rates in
        2040, 2050 and 2060.
        """
        names = {v: k for k, v in SCREENING_TESTS.items()}
        long, elim = [], []
        for unit in [INDIA, *STATES]:
            res = (self.national_test_exp if unit == INDIA else self.states_test_exp).result(unit)
            model = HPVCervicalCancerModel.for_unit(unit)
            cir0, cmr0 = rates(model.simulate(), model.params, model.inputs.denominator)
            years = model.years[1:].astype(int)
            base = {int(y): (cir0[i], cmr0[i]) for i, y in enumerate(years)}
            for si, sc in enumerate(res.scenarios):
                test = names[round(sc.H_pos, 2)]
                first, snap = None, {}
                for y in range(SCREENING_START, int(years[-1]) + 1):
                    s = res.snapshot(si, y)
                    b_cir, b_cmr = base[y]
                    if first is None and s.cir < WHO_THRESHOLD:
                        first = y
                    snap[y] = s
                    long.append(dict(
                        Unit=unit, Region=region(unit), Test_type=test, H_pos=sc.H_pos,
                        Coverage_pct=int(round(sc.coverage * 100)), Year=y,
                        CIR_mean=s.cir, CIR_lo95=max(s.cir - s.cir_halfwidth, 0), CIR_hi95=s.cir + s.cir_halfwidth,
                        CMR_mean=s.cmr, CMR_lo95=max(s.cmr - s.cmr_halfwidth, 0), CMR_hi95=s.cmr + s.cmr_halfwidth,
                        CIR_no_screening=b_cir, CMR_no_screening=b_cmr,
                        CIR_reduction_pct=100 * (1 - s.cir / b_cir) if b_cir > 0 else np.nan,
                        CMR_reduction_pct=100 * (1 - s.cmr / b_cmr) if b_cmr > 0 else np.nan))
                row = dict(Unit=unit, Region=region(unit), Test_type=test, H_pos=sc.H_pos,
                           Coverage_pct=int(round(sc.coverage * 100)), First_year_CIR_below4=first)
                for y in (2040, 2050, 2060):
                    s = snap[y]
                    row.update({f"CIR_{y}": s.cir, f"CIR_{y}_hw95": s.cir_halfwidth,
                                f"CMR_{y}": s.cmr, f"CMR_{y}_hw95": s.cmr_halfwidth,
                                f"CIR_no_screening_{y}": base[y][0], f"CMR_no_screening_{y}": base[y][1],
                                f"CIR_reduction_pct_{y}": 100 * (1 - s.cir / base[y][0]),
                                f"CMR_reduction_pct_{y}": 100 * (1 - s.cmr / base[y][1])})
                elim.append(row)
        el = pd.DataFrame(elim)
        el["First_year_CIR_below4"] = el["First_year_CIR_below4"].astype("Int64")
        return {"TestType_by_Unit": pd.DataFrame(long), "TestType_Elimination": el}

    def test_type_table(self) -> pd.DataFrame:
        """National test-type comparison (75 % coverage), annual rows 2003–2070 (Table S8).

        Intervals follow the convention of every other table and figure (mean ± half the
        spread of the 2.5/97.5-percentile trajectories); the no-screening reference is the
        deterministic run.
        """
        res = self.test_exp.result(INDIA)
        model = HPVCervicalCancerModel.for_unit(INDIA)
        p = model.params
        base = model.simulate()
        N, ys = res.denominator, res.years
        names = {v: k for k, v in SCREENING_TESTS.items()}

        def infl(cp, cc):
            return max((cc - cp) + (p.r_C + p.mu) * cp + p.mu_C * max(cp, 0) ** p.n_PL, 0.0)

        rows = []
        for si, sc in enumerate(res.scenarios):
            C, D = res.stats(si, "C"), res.stats(si, "C_D")
            for i in range(1, len(ys)):
                Nm = N[i]
                cm, cl, ch = infl(C.mean[i - 1], C.mean[i]), infl(C.low[i - 1], C.low[i]), infl(C.high[i - 1], C.high[i])
                cb = infl(base[i - 1, 6], base[i, 6])
                dm = max(D.mean[i] - D.mean[i - 1], 0)
                dl = max(D.low[i] - D.low[i - 1], 0)
                dh = max(D.high[i] - D.high[i - 1], 0)
                db = max(base[i, 9] - base[i - 1, 9], 0)
                rows.append(dict(
                    Scenario="Test-type sweep (coverage 75%)", Coverage_pct=int(round(sc.coverage * 100)),
                    H_pos=sc.H_pos, Test_type=names[round(sc.H_pos, 2)], Year=int(ys[i]),
                    CIR_mean=cm / Nm * 1e5, CIR_lo95=max(cm - abs(ch - cl) / 2, 0) / Nm * 1e5,
                    CIR_hi95=(cm + abs(ch - cl) / 2) / Nm * 1e5,
                    CMR_mean=dm / Nm * 1e5, CMR_lo95=max(dm - abs(dh - dl) / 2, 0) / Nm * 1e5,
                    CMR_hi95=(dm + abs(dh - dl) / 2) / Nm * 1e5,
                    Cases_mean=cm * 1000, Deaths_mean=dm * 1000,
                    Cases_averted_vs_noscreening=(cb - cm) * 1000, Deaths_averted_vs_noscreening=(db - dm) * 1000))
        return pd.DataFrame(rows)

    def save(self, frames: dict[str, pd.DataFrame] | None = None, directory: Path | None = None) -> Path:
        frames = frames or self.build()
        out = ensure(Path(directory) if directory else TABLES_DIR / "csv")
        for name, df in frames.items():
            df.to_csv(out / f"{name}.csv", index=False)
        return out
