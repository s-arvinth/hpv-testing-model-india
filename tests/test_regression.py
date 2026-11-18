"""Regenerated ensembles, tables and the high-level API vs the reference results in
paper/reference (slow: runs the paper's experiments unless they are cached)."""
import numpy as np
import pandas as pd
import pytest

from hpvcc.paper.experiments import Experiment
from hpvcc.paths import REFERENCE_DIR

REFERENCE = REFERENCE_DIR / "reference_numbers.csv"
pytestmark = pytest.mark.slow


@pytest.mark.parametrize("unit,experiment", [("Tamil Nadu", "states_coverage_sweep"),
                                             ("Bihar", "states_coverage_sweep"),
                                             ("India", "india_coverage_sweep")])
def test_coverage_sweep_matches_reference(unit, experiment):
    res = Experiment.get(experiment).result(unit)
    ref = pd.read_csv(REFERENCE)
    ref = ref[ref.State == unit]
    assert len(ref)
    for _, r in ref.iterrows():
        s = res.snapshot(res.index(coverage=r.Coverage_pct / 100), int(r.Year))
        assert abs(round(s.cir, 3) - r.CIR) <= 0.0005 + 1e-9
        assert abs(round(s.cmr, 3) - r.CMR) <= 0.0005 + 1e-9
        assert abs(round(s.cir_halfwidth, 3) - r.CIR_95CI_hw) <= 0.0005 + 1e-9


@pytest.mark.parametrize("unit", ["India", "Tamil Nadu"])
def test_project_reproduces_reference_numbers(unit):
    """The high-level API with its defaults gives the paper's numbers: draws are keyed by
    coverage, so any subset of coverages reproduces the full sweep."""
    from hpvcc.projection import project
    ref = pd.read_csv(REFERENCE)
    ref = ref[(ref.State == unit) & ref.Coverage_pct.isin([0, 25, 70])]
    proj = project(unit, coverage=[0.7, 0.25])
    for _, r in ref.iterrows():
        s = proj.snapshot(r.Coverage_pct / 100, int(r.Year))
        assert abs(round(s.cir, 3) - r.CIR) <= 0.0005 + 1e-9 and abs(round(s.cmr, 3) - r.CMR) <= 0.0005 + 1e-9
        assert abs(round(s.cir_halfwidth, 3) - r.CIR_95CI_hw) <= 0.0005 + 1e-9


def test_more_sensitive_test_gives_lower_incidence():
    res = Experiment.get("india_test_type").result("India")
    cir = [res.snapshot(i, 2050).cir for i in range(4)]        # Pap, VIA, LBC, HPV-DNA
    assert cir[0] > cir[1] > cir[2] > cir[3]


def test_latex_tables_match_reference(tmp_path):
    from hpvcc.paper.latex_tables import write_latex_tables
    from hpvcc.paper.tables import MasterTables
    frames = MasterTables().build()
    keys = {"A": "Annual_Long", "E": "Elimination_Year", "R": "Coverage_Required", "C": "Cumulative_Summary",
            "K": "Calibration_2012_2019", "P": "Parameters", "T": "India_TestType", "X": "TestType_Elimination"}
    paths = write_latex_tables({k: frames[v] for k, v in keys.items()}, tmp_path)
    for name, p in paths.items():
        assert p.read_text() == (REFERENCE_DIR / "paper_tables" / name).read_text(), name
