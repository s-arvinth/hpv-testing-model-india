"""Inflow incidence and the figure interval convention."""
import numpy as np

from hpvcc.metrics import TrajectoryStats, incidence_counts, snapshot
from hpvcc.parameters import load_parameters


def test_inflow_counts_cancers_that_leave_within_the_year():
    p = load_parameters("India")
    C = np.array([100.0, 100.0])               # stock unchanged ...
    out = (p.r_C + p.mu) * 100.0 + p.mu_C * 100.0 ** p.n_PL
    assert np.isclose(incidence_counts(C, p)[0], out)   # ... but cases = everything that left


def test_snapshot_halfwidth_is_half_the_percentile_spread():
    p = load_parameters("India")
    years = np.array([2049, 2050])
    C = TrajectoryStats(np.array([10.0, 10.5]), np.array([9.0, 9.2]), np.array([11.0, 11.9]))
    D = TrajectoryStats(np.array([0.0, 1.0]), np.array([0.0, 0.8]), np.array([0.0, 1.4]))
    s = snapshot(C, D, years, 2050, p, np.array([1e3, 1e3]))
    lo = incidence_counts(np.array([9.0, 9.2]), p)[0]
    hi = incidence_counts(np.array([11.0, 11.9]), p)[0]
    assert np.isclose(s.cir_halfwidth, abs(hi - lo) / 2 / 1e3 * 1e5)
    assert np.isclose(s.cmr_halfwidth, 0.3 / 1e3 * 1e5)
