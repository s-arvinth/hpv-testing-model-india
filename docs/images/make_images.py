"""Regenerate the figures shown in the README (India, paper defaults).

    python docs/images/make_images.py
"""
from pathlib import Path

from hpvcc import project
from hpvcc.viz.plots import plot_fit

HERE = Path(__file__).resolve().parent

proj = project("India", coverage=[0, 0.25, 0.5, 0.7, 0.9, 1.0], test="hpv-dna")
proj.plot_bars(HERE / "india_coverage_bars.png")
plot_fit("India", HERE / "india_fit.png")
print("README figures written to", HERE)
