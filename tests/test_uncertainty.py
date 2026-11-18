"""The shipped parameter sets regenerate the national uncertainty band (Fig S4)."""
import numpy as np

from hpvcc.model import HPVCervicalCancerModel
from hpvcc.paths import REFERENCE_DIR
from hpvcc.uncertainty import ParameterSetEnsemble


def test_parameter_sets_reproduce_monte_carlo_band():
    ref = np.load(REFERENCE_DIR / "national_fit_reference.npz")
    res = ParameterSetEnsemble.from_csv(HPVCervicalCancerModel.for_unit("India")).run()
    assert res.meta["n_ok"] == 500
    # tolerance covers platform numerics (~1e-7 relative)
    np.testing.assert_allclose(res.cir_low, ref["mc_CIR_l_mc"], atol=1e-4)
    np.testing.assert_allclose(res.cir_high, ref["mc_CIR_u_mc"], atol=1e-4)
    np.testing.assert_allclose(res.cir_best, ref["mc_CIR_best_mc"], atol=1e-4)
    np.testing.assert_allclose(res.cmr_low, ref["mc_CMR_l_mc"], atol=1e-4)
