"""Repository locations.

Everything is resolved relative to the repository root, so the code runs unchanged on any
machine. Set ``HPVCC_RESULTS`` to write outputs somewhere other than ``results/``.
"""
from __future__ import annotations

import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

# Model definition and inputs
PARAMS_DIR = ROOT / "params"                 #: fitted parameter sets (CSV)
DATA_RAW = ROOT / "data" / "raw"             #: population and prevalence sources
DATA_PROCESSED = ROOT / "data" / "processed" #: per-unit model inputs (JSON)
OBSERVED_RATES_FILE = ROOT / "data" / "registry" / "observed_rates_2012_2019.csv"  #: fit targets (shown, not fitted)

# Figure assets
ASSETS = ROOT / "assets"
FONTS_DIR = ASSETS / "fonts"
GEO_DIR = ASSETS / "geo"
DIAGRAMS_DIR = ASSETS / "diagrams"

# Reproduction of the paper (kept apart from the model itself)
PAPER_DIR = ROOT / "paper"
EXPERIMENTS_FILE = PAPER_DIR / "experiments.yaml"
REFERENCE_DIR = PAPER_DIR / "reference"      #: the paper's tables, figures and numbers

# Outputs
RESULTS = Path(os.environ.get("HPVCC_RESULTS", ROOT / "results"))
PAPER_RESULTS = RESULTS / "paper"
SIMULATIONS_DIR = PAPER_RESULTS / "simulations"   #: cached paper ensembles (.npz)
TABLES_DIR = PAPER_RESULTS / "tables"             #: master CSV tables + LaTeX tables
FIGURES_DIR = PAPER_RESULTS / "figures"           #: paper figures
FIGURES_TRANSPARENT_DIR = PAPER_RESULTS / "figures_transparent"  #: the same, without background


def ensure(path: Path) -> Path:
    """Create ``path`` (a directory) if needed and return it."""
    path.mkdir(parents=True, exist_ok=True)
    return path
