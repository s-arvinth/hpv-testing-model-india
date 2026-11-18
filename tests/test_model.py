"""Deterministic model: reference baseline values and basic invariants."""
from pathlib import Path

import numpy as np
import pytest

from hpvcc.metrics import rates
from hpvcc.model import HPVCervicalCancerModel
from hpvcc.paths import REFERENCE_DIR
from hpvcc.units import STATES


def _cir(unit, year, coverage=0.0):
    """Deterministic CIR; screened runs use the midpoints of the sampled rates
    (r_H = (r_I + 1) / 2, theta = gamma / 2)."""
    m = HPVCervicalCancerModel.for_unit(unit)
    traj = m.simulate(coverage, *m.midpoint_rates()) if coverage else m.simulate()
    cir, cmr = rates(traj, m.params, m.inputs.denominator)
    return cir[list(m.years[1:].astype(int)).index(year)]


def test_national_baseline_2050_is_11_56():
    """Headline no-screening CIR of the national model (paper: 11.56 per 100,000)."""
    assert round(_cir("India", 2050), 2) == 11.56


def test_national_baseline_matches_reference():
    ref = np.load(REFERENCE_DIR / "national_fit_reference.npz")["traj_baseline"]
    traj = HPVCervicalCancerModel.for_unit("India").simulate()
    np.testing.assert_allclose(traj, ref, rtol=1e-12, atol=0)


def test_screening_lowers_incidence():
    for unit in ("India", "Tamil Nadu", "Bihar"):
        assert _cir(unit, 2050, 0.7) < _cir(unit, 2050, 0.0)


def test_no_screening_leaves_detected_compartment_empty():
    traj = HPVCervicalCancerModel.for_unit("Punjab").simulate()
    assert np.all(traj[:, 4] == 0)


@pytest.mark.parametrize("unit", STATES[:6])
def test_states_rates_positive(unit):
    m = HPVCervicalCancerModel.for_unit(unit)
    cir, cmr = rates(m.simulate(), m.params, m.inputs.denominator)
    assert np.all(cir >= 0) and np.all(cmr >= 0)


def test_screened_run_needs_the_sampled_rates():
    """r_H and theta are sampled, not stored: a deterministic screened run must be given them."""
    m = HPVCervicalCancerModel.for_unit("Goa")
    assert not hasattr(m.params, "theta") and not hasattr(m.params, "r_H")
    with pytest.raises(ValueError):
        m.simulate(0.7)


def test_missed_infections_stay_infected():
    """False negatives join the undetected pool: with a useless test (H+ = 0) screening
    changes nothing, and a more sensitive test gives lower incidence."""
    m = HPVCervicalCancerModel.for_unit("Goa")
    rates_mid = m.midpoint_rates()
    blind = m.with_params(H_pos=0.0).simulate(0.9, *rates_mid)
    np.testing.assert_allclose(blind, m.simulate(), rtol=1e-6, atol=1e-9)
    c = [m.with_params(H_pos=h).simulate(0.75, *rates_mid)[-1, 9] for h in (0.70, 0.75, 0.85, 0.95)]
    assert c[0] > c[1] > c[2] > c[3]          # cumulative cancer deaths: Pap > VIA > LBC > HPV-DNA


def test_screening_conserves_people():
    """Screening only redirects infections; the living population is unchanged at t = start."""
    from hpvcc.model import rhs_no_screening, rhs_screening
    m = HPVCervicalCancerModel.for_unit("Goa")
    y = m.simulate()[30]
    a = np.array(rhs_no_screening(y, 2040.0, m._P, m.inputs.inflow))
    b = np.array(rhs_screening(y, 2040.0, m._P, m.inputs.inflow, 0.7, 2026, 0.6, 1e-4))
    assert abs(a[:3].sum() - b[:3].sum()) < 1e-9 and abs(a[3] - b[3]) < 1e-12   # same outflow from S
    assert abs((a[4] + a[5]) - (b[4] + b[5])) < 1e-9                              # same total new infections
