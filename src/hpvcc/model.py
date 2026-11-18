r"""
The 10-compartment HPV / cervical-cancer ODE system.

State vector (thousands of women)::

    S18  S30  S61   susceptible, ages 18–29 / 30–65 / 66–85
    S_R             re-susceptible after loss of immunity
    H_T             screen-detected HPV infection (tested positive)
    I               undetected HPV infection
    C               cervical cancer
    R_H             recovered from HPV
    C_S             cancer survivors
    C_D             cumulative cancer deaths (counter)

Force of infection λ = β (I + H_T) / N for the three age bands and λ_R = β_R (I + H_T) / N
for S_R, with N the living population. Screening at coverage α (from 2026) diverts a
fraction α·H⁺ of new infections among 30–65-year-olds (and the matching share of
re-infections) into H_T; every other new infection — including those the test misses —
enters I. Full equations: ``docs/model.md``.

The clearance rate r_H and the progression rate θ of screen-detected infection are not
parameters of a unit: they are sampled (r_H ~ U(r_I, 1), θ ~ U(0, γ), redrawn every year;
:mod:`hpvcc.simulation`) and passed to the screened right-hand side explicitly.

Implementation note
-------------------
The right-hand sides have a fixed floating-point evaluation order, so that fitted
trajectories are reproduced exactly. Re-run the regression tests (``pytest``) after
touching the expressions.
"""
from __future__ import annotations

import warnings

import numpy as np
from scipy.integrate import ODEintWarning, odeint

from .constants import SCREENING_START
from .inputs import RecruitmentInflow, UnitInputs
from .parameters import ModelParameters

# A few stochastic draws (fast clearance of screen-detected infection) make LSODA report
# "excess work" in single years; the step still completes.
warnings.filterwarnings("ignore", message="Excess work done", category=ODEintWarning)

# Order of the parameter tuple passed to the right-hand sides.
_FIELDS = ("p", "K1", "K2", "beta", "beta_reinfection", "gamma", "r_I",
           "r_C", "mu_C", "n_PL", "omega", "epsilon", "mu", "H_pos")


def pack(params: ModelParameters) -> tuple:
    """Parameters as a flat tuple (fast to unpack inside the ODE right-hand side)."""
    return tuple(float(getattr(params, f)) for f in _FIELDS)


def rhs_no_screening(y, t, P, inflow: RecruitmentInflow):
    """dy/dt without screening: every infection is undetected (I) and H_T stays empty, so
    the screened-pathway rates r_H and θ do not appear."""
    p, K1, K2, beta, betaR, gamma, r_I, r_C, mu_C, n_PL, omega, eps, mu, _ = P
    S18, S30, S61, S_R, H_T, I, C, R_H, C_S, C_D = y
    N = max(S18 + S30 + S61 + S_R + H_T + I + C + R_H + C_S, 1e-12)
    lam = beta * (I + H_T) / N
    lamR = betaR * (I + H_T) / N
    recruit = p * inflow(t)
    death = mu_C * (max(C, 0.0) ** n_PL)
    return [recruit - (lam + K1 + mu) * S18,
            K1 * S18 - (lam + K2 + mu) * S30,
            K2 * S30 - (lam + mu) * S61,
            eps * R_H - (lamR + mu) * S_R,
            -mu * H_T,
            lam * S18 + lam * S30 + lam * S61 + lamR * S_R - (gamma + r_I + mu) * I,
            gamma * I + omega * C_S - (r_C + mu) * C - death,
            r_I * I - (eps + mu) * R_H,
            r_C * C - (omega + mu) * C_S,
            death]


def rhs_screening(y, t, P, inflow: RecruitmentInflow, coverage: float,
                  screening_start: float, r_H: float, theta: float):
    """dy/dt with a screening programme of coverage ``coverage`` from ``screening_start``.

    A screened woman's new infection is detected with probability ``H`` (the test
    sensitivity) and enters H_T; an infection the test misses is a false negative and joins
    the undetected pool I, like the infections of women who are not screened. ``r_H`` and
    ``theta`` are sampled by the stochastic ensembles (see :mod:`hpvcc.simulation`) and
    passed in; the other rates come from ``P``.
    """
    p, K1, K2, beta, betaR, gamma, r_I, r_C, mu_C, n_PL, omega, eps, mu, H = P
    S18, S30, S61, S_R, H_T, I, C, R_H, C_S, C_D = y
    N = max(S18 + S30 + S61 + S_R + H_T + I + C + R_H + C_S, 1e-12)
    lam = beta * (I + H_T) / N
    lamR = betaR * (I + H_T) / N
    a = 0.0 if t < screening_start else coverage
    # Re-infections are screened in proportion to the 30–65 band's share of new infections.
    inflow_30 = (1 - a) * lam * S30
    total = inflow_30 + lam * S18 + lam * S61 + lamR * S_R
    frac_30 = inflow_30 / total if total > 0 else 0.0
    aR = a * frac_30
    recruit = p * inflow(t)
    death = mu_C * (max(C, 0.0) ** n_PL)
    return [recruit - (lam + K1 + mu) * S18,
            K1 * S18 - (lam + K2 + mu) * S30,
            K2 * S30 - (lam + mu) * S61,
            eps * R_H - (lamR + mu) * S_R,
            a * H * lam * S30 + aR * H * lamR * S_R - (theta + r_H + mu) * H_T,
            (1 - a * H) * lam * S30 + lam * S18 + lam * S61 + (1 - aR * H) * lamR * S_R - (gamma + r_I + mu) * I,
            theta * H_T + gamma * I + omega * C_S - (r_C + mu) * C - death,
            r_H * H_T + r_I * I - (eps + mu) * R_H,
            r_C * C - (omega + mu) * C_S,
            death]


class HPVCervicalCancerModel:
    """The ODE model for one unit (parameters + inputs).

    The system is integrated one calendar year at a time (``odeint`` restarted every
    year); this is what lets the stochastic
    ensembles redraw rates between years.

    Parameters
    ----------
    params : ModelParameters
        Rates of the unit (see :mod:`hpvcc.parameters`).
    inputs : UnitInputs
        Population, recruitment and initial conditions (see :mod:`hpvcc.inputs`).
    screening_start : float
        First year of the screening programme (default 2026).

    Example
    -------
    >>> model = HPVCervicalCancerModel.for_unit("India")
    >>> traj = model.simulate()                          # no screening, shape (69, 10)
    >>> faster = model.with_params(gamma=2 * model.params.gamma)   # what-if: faster progression
    >>> traj_fast = faster.simulate()

    Screened runs are normally done as stochastic ensembles (:mod:`hpvcc.simulation`), or
    through the high-level :func:`hpvcc.projection.project`.
    """

    def __init__(self, params: ModelParameters, inputs: UnitInputs,
                 screening_start: float = SCREENING_START):
        self.params = params
        self.inputs = inputs
        self.screening_start = float(screening_start)
        self._P = pack(params)

    @classmethod
    def for_unit(cls, unit: str, params: ModelParameters | None = None, **kwargs) -> "HPVCervicalCancerModel":
        """Model of one unit (``"India"`` or a state/UT) with its fitted parameters.

        ``params`` replaces the fitted parameters; ``screening_start`` may be passed too.
        """
        from .parameters import load_parameters
        from .units import resolve_unit
        unit = resolve_unit(unit)
        return cls(params or load_parameters(unit), UnitInputs.load(unit), **kwargs)

    @property
    def years(self) -> np.ndarray:
        return self.inputs.years

    @property
    def unit(self) -> str:
        return self.inputs.unit

    def with_params(self, params: ModelParameters | None = None, **changes) -> "HPVCervicalCancerModel":
        """Copy of the model with other parameters: a full set and/or individual changes,
        e.g. ``model.with_params(beta=0.4, r_I=0.5)``."""
        params = (params or self.params).replace(**changes) if changes else (params or self.params)
        return HPVCervicalCancerModel(params, self.inputs, self.screening_start)

    def with_screening_start(self, year: float) -> "HPVCervicalCancerModel":
        """Copy of the model whose screening programme starts in ``year``."""
        return HPVCervicalCancerModel(self.params, self.inputs, year)

    # ------------------------------------------------------------------ one-year steps
    def step(self, y, t0: float, t1: float) -> np.ndarray:
        """Advance ``y`` from ``t0`` to ``t1`` without screening."""
        return odeint(rhs_no_screening, y, [t0, t1], args=(self._P, self.inputs.inflow))[-1]

    def step_screened(self, y, t0: float, t1: float, coverage: float,
                      r_H: float, theta: float, P: tuple | None = None) -> np.ndarray:
        """Advance ``y`` from ``t0`` to ``t1`` with screening (rates ``r_H``/``theta`` given)."""
        P = self._P if P is None else P
        return odeint(rhs_screening, y, [t0, t1],
                      args=(P, self.inputs.inflow, coverage, self.screening_start, r_H, theta))[-1]

    # ------------------------------------------------------------------ full runs
    def midpoint_rates(self) -> tuple[float, float]:
        """Midpoints of the sampling ranges: ``r_H = (r_I + 1) / 2``, ``θ = γ / 2``."""
        return (self.params.r_I + 1.0) / 2.0, self.params.gamma / 2.0

    def simulate(self, coverage: float = 0.0, r_H: float | None = None,
                 theta: float | None = None) -> np.ndarray:
        """Deterministic trajectory, shape ``(n_years, 10)``.

        With ``coverage == 0`` the no-screening system is used. A screened run needs the
        two sampled rates of the screened pathway: pass ``r_H`` and ``theta`` (e.g.
        ``*model.midpoint_rates()``), or use a stochastic ensemble.
        """
        if coverage != 0.0 and (r_H is None or theta is None):
            raise ValueError("r_H and theta are sampled, not stored: pass both for a deterministic "
                             "screened run (e.g. model.simulate(0.7, *model.midpoint_rates())) "
                             "or run a StochasticEnsemble / project()")
        ys = self.years
        traj = np.empty((len(ys), 10))
        y = np.array(self.inputs.y0, dtype=float)
        traj[0] = y
        for i in range(1, len(ys)):
            if coverage == 0.0:
                y = self.step(y, ys[i - 1], ys[i])
            else:
                y = self.step_screened(y, ys[i - 1], ys[i], coverage, r_H, theta)
            traj[i] = y
        return traj

    def simulate_single_call(self, params: ModelParameters | None = None) -> np.ndarray:
        """No-screening trajectory from ONE ``odeint`` call over the whole grid.

        Used by the parameter-uncertainty analysis (:mod:`hpvcc.uncertainty`), which
        integrates each parameter set this way.
        """
        P = pack(params) if params is not None else self._P
        return odeint(rhs_no_screening, np.array(self.inputs.y0, dtype=float), self.years,
                      args=(P, self.inputs.inflow))
