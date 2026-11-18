"""Processed inputs are present and are exactly what the raw data produce."""
import numpy as np
import pytest

from hpvcc.data.national import build_national_inputs
from hpvcc.data.sources import RawData
from hpvcc.data.states import build_state_inputs
from hpvcc.inputs import UnitInputs
from hpvcc.units import INDIA, STATES


@pytest.fixture(scope="module")
def raw():
    return RawData()


def test_all_units_have_inputs():
    for unit in [INDIA, *STATES]:
        inp = UnitInputs.load(unit)
        assert inp.years[0] == 2002 and inp.years[-1] == 2070
        assert inp.y0.shape == (10,) and np.all(inp.y0 >= 0)
        assert np.all(inp.denominator > 0)


@pytest.mark.parametrize("unit", [INDIA, "Tamil Nadu", "Dadra and Nagar Haveli and Daman and Diu", "Ladakh"])
def test_inputs_rebuild_identically(unit, raw):
    stored = UnitInputs.load(unit)
    built = build_national_inputs(raw) if unit == INDIA else build_state_inputs(unit, raw)
    np.testing.assert_array_equal(built.y0, stored.y0)
    np.testing.assert_array_equal(built.denominator, stored.denominator)
    np.testing.assert_array_equal(built.inflow.values, stored.inflow.values)
    for band in ("S18", "S30", "S61"):
        np.testing.assert_array_equal(built.population[band], stored.population[band])


def test_inputs_hold_only_model_inputs():
    """Model inputs are population, recruitment and initial conditions only (the observed
    registry rates live apart, in data/registry/)."""
    import json
    from hpvcc.inputs import UNITS_DIR
    for path in UNITS_DIR.glob("*.json"):
        assert set(json.loads(path.read_text())) == {"unit", "years", "population", "denominator",
                                                     "inflow", "y0", "provenance"}, path.name


def test_observed_rates_cover_every_unit():
    from hpvcc.observed import ObservedRates, model_fit
    from hpvcc.units import all_units
    for unit in all_units():
        obs = ObservedRates.load(unit)
        assert list(obs.years.astype(int)) == list(range(2012, 2020))
        assert (obs.cir > 0).all() and (obs.cmr > 0).all()
    fit = model_fit("Tamil Nadu")
    assert len(fit) == 8 and abs(fit.Model_CIR.mean() / fit.Observed_CIR.mean() - 1) < 0.1
