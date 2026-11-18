"""User-facing API: parameters, projections, command line."""
import numpy as np
import pandas as pd
import pytest

from hpvcc import HPVCervicalCancerModel, load_parameters, project
from hpvcc.cli import main
from hpvcc.units import resolve_unit


def test_unit_names_are_forgiving():
    assert resolve_unit("tamil_nadu") == "Tamil Nadu"
    assert resolve_unit("andaman & nicobar") == "Andaman and Nicobar"
    with pytest.raises(KeyError):
        resolve_unit("Atlantis")


def test_parameter_changes_and_custom_file(tmp_path):
    p = load_parameters("Goa")
    m = HPVCervicalCancerModel.for_unit("Goa").with_params(beta=2 * p.beta)
    assert m.params.beta == 2 * p.beta and m.params.gamma == p.gamma
    f = tmp_path / "my_params.csv"
    pd.DataFrame([{**p.as_dict(), "beta": 0.123}]).to_csv(f, index=False)
    assert load_parameters("Goa", file=f) == p.replace(beta=0.123)


def test_parameter_sets_hold_no_sampled_rates():
    """theta and r_H are sampled, so no parameter set (file or object) carries them."""
    for unit in ("India", "Tamil Nadu"):
        d = load_parameters(unit).as_dict()
        assert "theta" not in d and "r_H" not in d and d["H_pos"] == 0.95
    for f in ("india_final.csv", "states_final.csv", "india_uncertainty_sets.csv"):
        from hpvcc.paths import PARAMS_DIR
        cols = pd.read_csv(PARAMS_DIR / f, nrows=1).columns
        assert "theta" not in cols and "r_H" not in cols


def test_deterministic_projection_table():
    proj = project("Goa", coverage=[0, 0.7], deterministic=True)
    t = proj.table
    assert set(t.Coverage) == {0.0, 0.7} and t.Year.min() == 2003 and t.Year.max() == 2070
    assert (t[t.Coverage == 0].Test == "No screening").all()
    base, scr = (t[(t.Coverage == c) & (t.Year == 2050)].iloc[0] for c in (0.0, 0.7))
    assert scr.CIR < base.CIR
    assert list(proj.summary().columns[:2]) == ["Coverage", "Test"]


def test_small_stochastic_projection_has_intervals():
    proj = project("Goa", coverage=0.5, test=["via", "hpv-dna"], n_sims=20, workers=1)
    t = proj.table
    assert set(t.Test) == {"VIA", "HPV DNA test"}
    assert (t.CIR_hi95 >= t.CIR).all() and (t.CIR_lo95 <= t.CIR).all()
    assert proj.ensemble.n_sims == 20


def test_cli_run_and_params(tmp_path, capsys):
    out = tmp_path / "goa.csv"
    main(["run", "Goa", "--coverage", "0", "0.7", "--deterministic", "--out", str(out)])
    assert "CIR_2050" in capsys.readouterr().out
    assert len(pd.read_csv(out)) == 2 * 68
    main(["params", "Goa"])
    assert "beta" in capsys.readouterr().out


def test_plots_fit_and_background(tmp_path, capsys):
    import matplotlib.image as mpimg
    proj = project("Tamil Nadu", coverage=[0, 0.7], deterministic=True)
    assert proj.plot(tmp_path / "ts.png").exists()
    assert proj.plot_bars(tmp_path / "bars.png").exists()
    snap = proj.snapshot(0.7, 2050)
    row = proj.at(2050).iloc[1]
    assert abs(snap.cir - row.CIR) < 1e-9 and snap.cir_halfwidth == 0
    from hpvcc.viz.plots import plot_compartments, plot_fit
    assert plot_fit("Tamil Nadu", tmp_path / "fit.png").exists()
    assert plot_compartments(proj.model, path=tmp_path / "comp.png").exists()
    white = mpimg.imread(tmp_path / "fit.png")
    clear = mpimg.imread(plot_fit("Tamil Nadu", tmp_path / "fit_t.png", transparent=True))
    assert white.shape[-1] == 3 or white[0, 0, -1] == 1.0
    assert clear[0, 0, -1] == 0.0
    main(["fit", "Tamil Nadu"])
    assert "Observed_CIR" in capsys.readouterr().out
