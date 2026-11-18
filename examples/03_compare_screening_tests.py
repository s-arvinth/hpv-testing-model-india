"""
03 — Four screening tests at the same coverage (national model).

Tests differ only by their sensitivity H+ (Pap 0.70, VIA 0.75, LBC 0.85, HPV-DNA 0.95).
An infection the test misses stays undetected, so a more sensitive test gives lower
incidence. The second figure rebuilds the paper's Fig S3 (screen-detected infections over
time, CIR / CMR at three years); it needs every compartment, hence its own ensemble.

    python examples/03_compare_screening_tests.py
"""
from pathlib import Path

from hpvcc import HPVCervicalCancerModel, Scenario, project
from hpvcc.simulation import StochasticEnsemble, default_draws
from hpvcc.viz.figures import test_bars

OUT = Path("results/examples")
COVERAGE, TESTS = 0.75, ["pap", "via", "lbc", "hpv-dna"]

proj = project("India", coverage=COVERAGE, test=TESTS)
print(proj.summary(years=(2030, 2050, 2070)).round(2).to_string(index=False))
print(f"\nfigure -> {proj.plot(OUT / '03_india_tests.png')}")

# Fig S3: the same 1000 simulations per test, with every compartment kept (for H_T)
model = HPVCervicalCancerModel.for_unit("India")
scenarios = [Scenario.for_test(COVERAGE, t) for t in TESTS]
ens = StochasticEnsemble(model, scenarios, n_sims=1000, draws=default_draws(scenarios=scenarios),
                         keep_full=True).run()
years = (2030, 2050, 2070)
detected = [ens.stats(i, "H_T") for i in range(len(TESTS))]
snaps = [{Y: ens.snapshot(i, Y) for Y in years} for i in range(len(TESTS))]
fig = test_bars(ens.years, detected, snaps, [s.H_pos for s in ens.scenarios], OUT / "03_india_fig_s3.png", years)
print(f"figure -> {fig}")

# A sensitivity that is not one of the named tests
custom = project("India", coverage=COVERAGE, H_pos=0.90)
print(custom.summary(years=years).round(2).to_string(index=False))
