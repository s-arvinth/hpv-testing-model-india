"""
National-model figures of the appendix.

* :func:`test_type_figure` — Fig S3: screen-detected infections H_T over time (left) and
  CIR / CMR in 2030, 2050 and 2070 (right) for four screening tests at 75 % coverage.
* :func:`monte_carlo_bars` — Fig S4: best-fit no-screening CIR (or CMR) every ten years
  with the 95 % band of the accepted Monte Carlo (GLUE) parameter sets.
* :func:`model_diagram` — Figs S1, S2.
"""
from __future__ import annotations

from pathlib import Path

from ...units import INDIA
from ...viz.figures import TEST_LABELS, test_bars, uncertainty_bars  # noqa: F401
from ..experiments import Experiment


def test_type_figure(out: Path, years=(2030, 2050, 2070)) -> Path:
    """Fig S3 from the test-type ensemble (see :func:`hpvcc.viz.figures.test_bars`)."""
    res = Experiment.get("india_test_type").result(INDIA)
    sens = [s.H_pos for s in res.scenarios]
    detected = [res.stats(si, "H_T") for si in range(len(sens))]
    snaps = [{Y: res.snapshot(si, Y) for Y in years} for si in range(len(sens))]
    return test_bars(res.years, detected, snaps, sens, out, years)


def monte_carlo_bars(metric: str, out: Path, step: int = 10, width: float = 6.0) -> Path:
    """Fig S4 panel from the parameter-set ensemble."""
    g = Experiment.get("india_monte_carlo").result(INDIA)
    return uncertainty_bars(g, metric, out, step, width, n_label=g.meta.get("n_samples_drawn", len(g.sets)))


def model_diagram(name: str, out: Path, variant: str | None = None) -> Path:
    """Figs S1/S2: copy the model diagram (``name`` = "full" or "baseline"; ``variant`` =
    "paper", or "dark" for the white-line version on a transparent background; by default
    "dark" when figures are saved without a background)."""
    import shutil
    from ...paths import DIAGRAMS_DIR
    from ...viz.style import is_transparent
    variant = variant or ("dark" if is_transparent() else "paper")
    files = {"full": "HPV_Testing_Model_for_Paper", "baseline": "HPV_Testing_Baseline_Model_for_Paper"}
    src = DIAGRAMS_DIR / variant / (files[name] + ("_for_Git" if variant == "dark" else "") + ".pdf")
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(src, out)
    return out


__all__ = ["test_type_figure", "monte_carlo_bars", "model_diagram"]
