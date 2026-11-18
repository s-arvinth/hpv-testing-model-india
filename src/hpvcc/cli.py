"""
hpvcc — run the HPV / cervical-cancer model from the command line.

Model
    hpvcc units                                  list the national and 36 state/UT models
    hpvcc params India                           fitted parameters of one unit
    hpvcc fit India --plot fit.png               fitted model vs observed rates, 2012–2019
    hpvcc run India                              no-screening projection, 2003–2070
    hpvcc run India --coverage 0 0.5 0.7         HPV-DNA testing at 50 % and 70 % coverage from 2026
    hpvcc run Kerala --coverage 0 0.7            the same for a state/UT
    hpvcc run India --coverage 0.75 --test pap via lbc hpv-dna --sims 200
    hpvcc run Bihar --coverage 0.7 --set beta=0.5 gamma=0.001 --plot bihar.png --out bihar.csv
    hpvcc prepare-inputs                         rebuild data/processed from data/raw

Figures are drawn in the paper's style; `hpvcc --no-background ...` saves them without a
background.

Reproducing the paper (see paper/README.md)
    hpvcc paper list | simulate [EXPERIMENT ...] | tables | figure fig3 [figS11 ...] | all
"""
from __future__ import annotations

import argparse
import sys
import time
import warnings
from pathlib import Path


# ====================================================================== model commands
def _cmd_units(args):
    from .units import INDIA, REGIONS
    print(f"{INDIA}  (national model)")
    for region, members in REGIONS.items():
        print(f"\n{region}")
        for u in members:
            print(f"  {u}")


def _cmd_params(args):
    from .parameters import PARAMETER_INFO, SAMPLED_RATES, load_parameters
    p = load_parameters(args.unit, file=args.params_file)
    print(f"Parameters of {args.unit}" + (f" ({args.params_file})" if args.params_file else ""))
    for name, value in p.as_dict().items():
        sym, meaning, unit = PARAMETER_INFO[name]
        print(f"  {name:17s} {value:<14.6g} {sym:5s} {meaning} [{unit}]")
    for name, (sym, meaning, dist) in SAMPLED_RATES.items():
        print(f"  {name:17s} {'sampled':<14s} {sym:5s} {meaning} [{dist}]")


def _parse_changes(items):
    changes = {}
    for item in items or []:
        if "=" not in item:
            raise SystemExit(f"--set expects NAME=VALUE, got {item!r}")
        k, v = item.split("=", 1)
        changes[k.strip()] = float(v)
    return changes


def _cmd_run(args):
    import pandas as pd
    from .projection import project
    t0 = time.time()
    proj = project(args.unit, coverage=args.coverage, test=args.test, H_pos=args.sensitivity,
                   n_sims=args.sims, deterministic=args.deterministic, params_file=args.params_file,
                   changes=_parse_changes(args.set), screening_start=args.start, seed=args.seed,
                   workers=args.workers)
    with pd.option_context("display.width", 200, "display.max_columns", 50):
        print(f"\n{proj.unit}: mean CIR and CMR per 100,000 women "
              f"({', '.join(sorted(set(proj.table.Method)))})\n")
        print(proj.summary(args.years).round(2).to_string(index=False))
        if any(s.coverage > 0 for s in proj.scenarios):
            print("\nFirst year with mean CIR < 4 per 100,000:")
            print(proj.elimination_year().to_string(index=False))
    if args.out:
        print(f"\nannual table -> {proj.to_csv(args.out)}")
    if args.plot:
        print(f"figure -> {proj.plot(args.plot)}")
    if args.bars:
        print(f"figure -> {proj.plot_bars(args.bars, years=args.bar_years)}")
    print(f"({time.time() - t0:.0f} s)")


def _cmd_fit(args):
    import pandas as pd
    from .model import HPVCervicalCancerModel
    from .observed import model_fit, observed_label
    from .parameters import load_parameters
    from .units import resolve_unit
    unit = resolve_unit(args.unit)
    model = None
    if args.params_file or args.set:
        params = load_parameters(unit, file=args.params_file) if args.params_file else None
        model = HPVCervicalCancerModel.for_unit(unit, params=params)
        if args.set:
            model = model.with_params(**_parse_changes(args.set))
    f = model_fit(unit, model)
    with pd.option_context("display.width", 200):
        print(f"\n{unit}: fitted model (no screening) vs observed, {observed_label(unit)}, per 100,000 women\n")
        print(f.drop(columns="Unit").round(3).to_string(index=False))
        m = f.mean(numeric_only=True)
        print(f"\nmean 2012–2019: CIR {m.Model_CIR:.2f} model vs {m.Observed_CIR:.2f} observed "
              f"({100 * (m.Model_CIR / m.Observed_CIR - 1):+.1f} %), CMR {m.Model_CMR:.2f} vs "
              f"{m.Observed_CMR:.2f} ({100 * (m.Model_CMR / m.Observed_CMR - 1):+.1f} %)")
    if args.plot:
        from .viz.plots import plot_fit
        print(f"figure -> {plot_fit(unit, args.plot, model=model)}")


def _cmd_prepare(args):
    from .data import build_all
    paths = build_all(verbose=not args.quiet)
    print(f"wrote {len(paths)} unit input files")


# ====================================================================== paper commands
_FRAME_KEYS = {"A": "Annual_Long", "E": "Elimination_Year", "R": "Coverage_Required", "C": "Cumulative_Summary",
               "K": "Calibration_2012_2019", "P": "Parameters", "T": "India_TestType", "X": "TestType_Elimination"}


def _cmd_paper_simulate(args):
    from .paper.experiments import Experiment
    names = args.experiments or [e.name for e in Experiment.all()]
    for name in names:
        exp = Experiment.get(name)
        print(f"{name}: {exp.description}")
        for unit in (args.unit or exp.units):
            exp.run(unit, n_sims=args.n_sims, workers=args.workers)


def _cmd_paper_tables(args):
    from .paper.latex_tables import write_latex_tables
    from .paper.tables import MasterTables
    mt = MasterTables()
    frames = mt.build(verbose=not args.quiet)
    csv_dir = mt.save(frames, Path(args.out) / "csv" if args.out else None)
    tex = write_latex_tables({k: frames[v] for k, v in _FRAME_KEYS.items()},
                             Path(args.out) / "latex" if args.out else None)
    print(f"master tables -> {csv_dir}")
    for name, p in tex.items():
        print(f"  {name} -> {p}")


def _cmd_paper_figure(args):
    from .paper.figures.registry import make_figure
    for key in args.figures:
        make_figure(key, Path(args.out) if args.out else None, transparent=args.no_background)


def _cmd_paper_all(args):
    t0 = time.time()
    args.experiments, args.unit = None, None
    _cmd_paper_simulate(args)
    _cmd_paper_tables(args)
    from .paper.figures.registry import make_all
    make_all(Path(args.out) / "figures" if args.out else None, transparent=args.no_background)
    print(f"done in {(time.time() - t0) / 60:.1f} min")


def _cmd_paper_list(args):
    from .paper.experiments import Experiment
    from .paper.figures.registry import REGISTRY
    print("Experiments (paper/experiments.yaml)")
    for e in Experiment.all():
        print(f"  {e.name:24s} {e.description}")
    print("\nFigures")
    for f in REGISTRY.values():
        src = ", ".join(f.experiments) or "no simulation"
        print(f"  {f.key:7s} {f.number:7s} {len(f.panels):3d} panel(s)  [{src}]  {f.caption}")


# ====================================================================== parser
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="hpvcc", description="HPV infection, testing and cervical-cancer model "
                                "for India and its states/UTs.", epilog="Run `hpvcc <command> -h` for details.")
    p.add_argument("--workers", type=int, default=None, help="parallel processes (default: all cores but two)")
    p.add_argument("--quiet", action="store_true")
    p.add_argument("--no-background", action="store_true",
                   help="save figures without a background (paper figures go to results/paper/figures_transparent/)")
    sub = p.add_subparsers(dest="cmd", required=True, metavar="command")

    sub.add_parser("units", help="list the modelled units").set_defaults(func=_cmd_units)

    pa = sub.add_parser("params", help="show the fitted parameters of a unit")
    pa.add_argument("unit")
    pa.add_argument("--params-file", default=None, help="read parameters from this CSV instead")
    pa.set_defaults(func=_cmd_params)

    r = sub.add_parser("run", help="project CIR and CMR for a unit under screening scenarios",
                       formatter_class=argparse.RawDescriptionHelpFormatter,
                       description="Project crude cervical-cancer incidence (CIR) and mortality (CMR) per "
                                   "100,000 women, 2003–2070, from the fitted parameters.")
    r.add_argument("unit", help='"India" or a state/UT (see `hpvcc units`)')
    r.add_argument("--coverage", type=float, nargs="+", default=[0.0],
                   help="screening coverage(s) in [0, 1]; 0 = no screening (default)")
    g = r.add_mutually_exclusive_group()
    g.add_argument("--test", nargs="+", default=None, help="screening test(s): hpv-dna (default), lbc, via, pap")
    g.add_argument("--sensitivity", type=float, nargs="+", default=None, help="test sensitivity H+ instead of --test")
    r.add_argument("--start", type=float, default=2026, help="first year of screening (default 2026)")
    r.add_argument("--sims", type=int, default=None,
                   help="stochastic simulations per screened scenario (default 1000, as in the paper)")
    r.add_argument("--deterministic", action="store_true", help="single run at the midpoints of the r_H and theta ranges (fast)")
    r.add_argument("--seed", type=int, default=0, help="seed of the random draws (default 0, as in the paper)")
    r.add_argument("--set", nargs="+", metavar="NAME=VALUE", help="change parameters, e.g. --set beta=0.4 r_I=0.5")
    r.add_argument("--params-file", default=None, help="CSV with your own parameters (format of params/*.csv)")
    r.add_argument("--years", type=int, nargs="+", default=[2030, 2040, 2050, 2060, 2070],
                   help="years shown in the printed summary")
    r.add_argument("--out", default=None, help="write the annual table (CSV)")
    r.add_argument("--plot", default=None, help="save a figure of CIR and CMR over time (.png/.pdf/.svg)")
    r.add_argument("--bars", default=None,
                   help="save the paper's coverage figure: CIR and CMR by coverage at --bar-years (one test)")
    r.add_argument("--bar-years", type=int, nargs="+", default=[2040, 2050, 2060], help="years of --bars")
    r.set_defaults(func=_cmd_run)

    fi = sub.add_parser("fit", help="fitted model vs observed registry rates, 2012–2019")
    fi.add_argument("unit")
    fi.add_argument("--plot", default=None, help="save the paper's model-fit figure (.png/.pdf/.svg)")
    fi.add_argument("--set", nargs="+", metavar="NAME=VALUE", help="change parameters first")
    fi.add_argument("--params-file", default=None, help="CSV with your own parameters")
    fi.set_defaults(func=_cmd_fit)

    sub.add_parser("prepare-inputs", help="rebuild data/processed from data/raw").set_defaults(func=_cmd_prepare)

    pp = sub.add_parser("paper", help="reproduce the paper's simulations, tables and figures")
    psub = pp.add_subparsers(dest="paper_cmd", required=True, metavar="paper-command")
    psub.add_parser("list", help="experiments and figures").set_defaults(func=_cmd_paper_list)
    s = psub.add_parser("simulate", help="run experiments (all by default) and cache them")
    s.add_argument("experiments", nargs="*")
    s.add_argument("--unit", action="append", help="only this unit (repeatable)")
    s.add_argument("--n-sims", type=int, default=None, help="fewer simulations for a quick look")
    s.set_defaults(func=_cmd_paper_simulate)
    t = psub.add_parser("tables", help="master CSV tables and the LaTeX tables")
    t.add_argument("--out", default=None)
    t.set_defaults(func=_cmd_paper_tables)
    f = psub.add_parser("figure", help="draw figures, e.g. `hpvcc paper figure fig3 figS4`")
    f.add_argument("figures", nargs="+")
    f.add_argument("--out", default=None)
    f.set_defaults(func=_cmd_paper_figure)
    a = psub.add_parser("all", help="every simulation, table and figure")
    a.add_argument("--out", default=None)
    a.add_argument("--n-sims", type=int, default=None)
    a.set_defaults(func=_cmd_paper_all)
    return p


def main(argv=None):
    warnings.filterwarnings("ignore", message="Excess work done")
    args = build_parser().parse_args(argv)
    if args.workers:
        import os
        os.environ["HPVCC_WORKERS"] = str(args.workers)
    if args.no_background:
        from .viz.style import set_transparent
        set_transparent(True)
    args.func(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
