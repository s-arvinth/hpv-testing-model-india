"""
02 — HPV-DNA testing at increasing coverage (national model).

Each screened scenario is a stochastic ensemble: the clearance (r_H) and progression (θ)
rates of screen-detected infections are redrawn every year, giving a 95 % interval. With
the defaults (1000 simulations, seed 0) the numbers — and the bar figure — are
those of the paper (Fig 3). Takes about 20 s on 40 cores; pass n_sims=200 for a quick look.

    python examples/02_screening_coverage.py
"""
from pathlib import Path

from hpvcc import project

OUT = Path("results/examples")

proj = project("India", coverage=[0, 0.25, 0.5, 0.7, 0.9, 1.0], test="hpv-dna")

print(proj.summary().round(2).to_string(index=False))
print("\nFirst year with mean CIR below 4 per 100,000 (elimination proxy):")
print(proj.elimination_year().to_string(index=False))

proj.to_csv(OUT / "02_india_coverage.csv")
print(f"\nfigure -> {proj.plot(OUT / '02_india_projection.png')}")         # rates over time
print(f"figure -> {proj.plot_bars(OUT / '02_india_coverage_bars.png')}")   # = Fig 3 of the paper
