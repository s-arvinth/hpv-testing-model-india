"""
LaTeX tables of the paper (Table 1 and Appendix Tables S4–S11), generated from the master
tables (:mod:`hpvcc.paper.tables`). Each table is a stand-alone ``.tex`` fragment that the
manuscript pulls in with ``\\input{}``; nothing is edited by hand.

==============================  ===================  ===========================
file                            paper                label
==============================  ===================  ===========================
``T2_state_policy_coverage``    Table 1              ``tab:state_policy``
``T1_national_coverage``        Table S4             ``tab:national_coverage``
``T3_elimination_timeline``     Table S5             ``tab:elimination``
``S1_regional_summary``         Table S6             ``tab:app_regions``
``S2_incremental_benefit``      Table S7             ``tab:app_incremental``
``S3_modality_comparison``      Table S8             ``tab:app_modality``
``S4_calibration_fit``          Table S9             ``tab:app_calibration``
``S5_state_policy_hpv_dna``     Table S10            ``tab:app_state_policy_hpv``
``S6_state_tests_other``        Table S11            ``tab:app_state_tests_other``
==============================  ===================  ===========================

Conventions: India = the dedicated national model; 95 % CI = ensemble-mean trajectory ±
the half-width drawn as error bars in the figures; elimination proxy = mean CIR < 4 per
100,000 (the WHO threshold is defined on age-standardised incidence, which this model
cannot compute). ``clean=False`` adds a note to each caption naming its source tables.
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd

from ..paths import TABLES_DIR, ensure

DECIMAL = '.'          # '.' for now; set to '\\textperiodcentered' / '·' for Lancet style
CLEAN = True           # False: add a source note to each caption
NEVER = '--'            # not reached by 2070
COVS_MAIN = [0, 25, 50, 70, 90, 100]
SWEEP_YEARS = [2040, 2050, 2060]   # snapshot years of the coverage-sweep figures
LATER = SWEEP_YEARS[1:]            # years shown for reductions / no-screening / CMR columns
YEARS = SWEEP_YEARS                # national, state and regional tables follow the sweep figures
MODALITY_YEARS = [2030, 2050, 2070]  # the screening-test figure (H+ sweep) uses these years
def yh(ys):
    return ' & '.join(str(y) for y in ys)
POLICY_COV = 70
REGION_ORDER = ['Central', 'East', 'North', 'Northeast', 'South', 'West']
INDIA = 'India'

# ----------------------------------------------------------------------------- formatting helpers
def fnum(x, nd=1):
    if x is None or (isinstance(x, float) and (math.isnan(x) or math.isinf(x))):
        return NEVER
    s = f"{x:,.{nd}f}"
    s = s.replace(',', '{,}')
    return s.replace('.', DECIMAL) if DECIMAL != '.' else s

def fint(x):
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return NEVER
    return f"{int(round(x)):,}".replace(',', '{,}')

def fk(x):
    """thousands, no decimals"""
    return NEVER if (x is None or (isinstance(x, float) and math.isnan(x))) else f"{x/1000:,.0f}".replace(',', '{,}')

def fci(m, lo, hi, nd=1):
    """mean (lo--hi); collapses to the mean when the interval is degenerate (deterministic arm)."""
    if fnum(lo, nd) == fnum(hi, nd) == fnum(m, nd):
        return fnum(m, nd)
    return f"{fnum(m, nd)} ({fnum(lo, nd)}--{fnum(hi, nd)})"

def fpct(x, nd=0):
    return NEVER if (x is None or (isinstance(x, float) and math.isnan(x))) else fnum(x, nd)

def fyear(y):
    return NEVER if (y is None or (isinstance(y, float) and math.isnan(y))) else str(int(y))

def tex_escape(s):
    return str(s).replace('&', '\\&').replace('%', '\\%')

# ----------------------------------------------------------------------------- LaTeX scaffolding
TRACK_OPEN = ""

def table_env(body, caption, label, landscape=False, size='\\footnotesize', note=None, width=None):
    env = 'sidewaystable' if landscape else 'table'
    cap = caption
    if note and not CLEAN:
        cap += "\\cnote{" + note + "}"
    out = [f"\\begin{{{env}}}[p]" if landscape else f"\\begin{{{env}}}[!htbp]",
           "\\centering", TRACK_OPEN.rstrip('\n'), size,
           f"\\caption{{{cap}}}", f"\\label{{{label}}}", body, f"\\end{{{env}}}"]
    return "\n".join(x for x in out if x) + "\n"

def tabular(colspec, header_rows, body_rows, footer=None):
    lines = [f"\\begin{{tabular}}{{{colspec}}}", "\\toprule"]
    lines += header_rows
    lines.append("\\midrule")
    lines += body_rows
    lines.append("\\bottomrule")
    lines.append("\\end{tabular}")
    return "\\begin{adjustbox}{max width=\\linewidth}\n" + "\n".join(lines) + "\n\\end{adjustbox}"

# ----------------------------------------------------------------------------- lookups
def row(A, unit, cov, year):
    r = A[(A.Unit == unit) & (A.Coverage_pct == cov) & (A.Year == year)]
    assert len(r) == 1, (unit, cov, year, len(r))
    return r.iloc[0]

def state_units(A):
    return sorted(u for u in A.Unit.unique() if u != INDIA)

def region_of(A):
    return A.drop_duplicates('Unit').set_index('Unit').Region.to_dict()

# ============================================================================= T1 national coverage x horizon
def t1_national(A, C, E):
    hdr = ["Coverage (\\%) & \\multicolumn{3}{c}{Crude incidence rate (95\\% CI)} & \\multicolumn{2}{c}{CIR reduction (\\%)} & "
           "\\multicolumn{3}{c}{Crude mortality rate (95\\% CI)} & \\multicolumn{2}{c}{CMR reduction (\\%)} & "
           "First year & \\multicolumn{4}{c}{Averted 2026--2070, thousands (95\\% CI)} \\\\",
           "\\cmidrule(lr){2-4}\\cmidrule(lr){5-6}\\cmidrule(lr){7-9}\\cmidrule(lr){10-11}\\cmidrule(lr){13-16}",
           f" & {yh(YEARS)} & {yh(LATER)} & {yh(YEARS)} & {yh(LATER)} & CIR $<$ 4 & Cases & \\% & Deaths & \\% \\\\"]
    body = []
    for cov in COVS_MAIN:
        rows = {y: row(A, INDIA, cov, y) for y in YEARS}
        cells = [f"{cov}"]
        cells += [fci(rows[y].CIR_mean, rows[y].CIR_lo95, rows[y].CIR_hi95) for y in YEARS]
        cells += [fpct(rows[y].CIR_reduction_pct_vs_baseline) for y in LATER]
        cells += [fci(rows[y].CMR_mean, rows[y].CMR_lo95, rows[y].CMR_hi95) for y in YEARS]
        cells += [fpct(rows[y].CMR_reduction_pct_vs_baseline) for y in LATER]
        e = E[(E.Unit == INDIA) & (E.Coverage_pct == cov)].iloc[0]
        cells.append(fyear(e['First_year_CIR_below4']))
        c = C[(C.Unit == INDIA) & (C.Coverage_pct == cov) & (C['Horizon_year(2026..Y)'] == 2070)].iloc[0]
        if cov == 0:
            cells += ['0', '0', '0', '0']
        else:
            cells.append(f"{fk(c.CumCases_averted_mean)} ({fk(c.CumCases_averted_lo95)}--{fk(c.CumCases_averted_hi95)})")
            cells.append(fpct(c.CumCases_averted_pct))
            cells.append(f"{fk(c.CumDeaths_averted_mean)} ({fk(c.CumDeaths_averted_lo95)}--{fk(c.CumDeaths_averted_hi95)})")
            cells.append(fpct(c.CumDeaths_averted_pct))
        body.append(" & ".join(cells) + " \\\\")
    base2026 = row(A, INDIA, 0, 2026)
    tab = tabular("l" + "c" * 15, hdr, body)
    caption = ("Projected cervical cancer incidence and mortality in India under a single-round HPV DNA testing programme "
               "for women aged 30--65 years, by programme coverage, at three snapshot years. Rates are per 100{,}000 women; "
               "reductions are relative to no screening in the same year; cases and deaths averted (in thousands) are cumulative from the start "
               "of screening in 2026 to 2070 (percentage of the no-screening cumulative total in parentheses). "
               f"No-screening CIR and CMR in 2026: {fnum(base2026.CIR_mean)} and {fnum(base2026.CMR_mean)} per 100{{,}}000. "
               "95\\% CI from 1000 stochastic simulations in which $r_H$ and $\\theta$ are redrawn every year (dedicated national model); "
               "the no-screening arm is deterministic. CIR=crude incidence rate. CMR=crude mortality rate. "
               "First year CIR $<$ 4 is the first calendar year in which the mean CIR falls below 4 per 100{,}000, used as a proxy for the WHO elimination threshold, which is defined on age-standardised incidence.")
    note = ("Generated by hpvcc (paper tables) from Annual\\_Long, Elimination\\_Year and Cumulative\\_Summary; dedicated national model "
            "(no-screening CIR 11.56 in 2050).")
    return table_env(tab, caption, "tab:national_coverage", landscape=True, size='\\footnotesize', note=note)

# ============================================================================= S5 state-by-state, HPV DNA testing, three years
def s5_states_hpv_dna(A, E, R, K):
    reg = region_of(A)
    units = state_units(A)
    # observed mean CIR/CMR 2012-2019
    kobs = K[K.Unit != INDIA].groupby('Unit')[['Observed_CIR', 'Observed_CMR']].mean()
    hdr = ["State / Union Territory & \\multicolumn{2}{c}{Observed 2012--19} & \\multicolumn{2}{c}{CIR, no screening} & "
           f"\\multicolumn{{3}}{{c}}{{CIR, {POLICY_COV}\\% coverage (95\\% CI)}} & \\multicolumn{{2}}{{c}}{{Reduction (\\%)}} & "
           f"\\multicolumn{{2}}{{c}}{{CMR, {POLICY_COV}\\% coverage}} & First year \\\\",
           "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\\cmidrule(lr){6-8}\\cmidrule(lr){9-10}\\cmidrule(lr){11-12}",
           f" & CIR & CMR & {yh(LATER)} & {yh(YEARS)} & {yh(LATER)} & {yh(LATER)} & CIR $<$ 4 \\\\"]
    body = []
    def unit_row(u, bold=False):
        cells = [("\\textbf{" + tex_escape(u) + "}") if bold else tex_escape(u)]
        if u in kobs.index:
            cells += [fnum(kobs.loc[u, 'Observed_CIR']), fnum(kobs.loc[u, 'Observed_CMR'])]
        else:
            ki = K[K.Unit == u][['Observed_CIR', 'Observed_CMR']].mean()
            cells += [fnum(ki.Observed_CIR), fnum(ki.Observed_CMR)]
        cells += [fnum(row(A, u, 0, y).CIR_mean) for y in LATER]
        rr = {y: row(A, u, POLICY_COV, y) for y in YEARS}
        cells += [fci(rr[y].CIR_mean, rr[y].CIR_lo95, rr[y].CIR_hi95) for y in YEARS]
        cells += [fpct(rr[y].CIR_reduction_pct_vs_baseline) for y in LATER]
        cells += [fnum(rr[y].CMR_mean) for y in LATER]
        e = E[(E.Unit == u) & (E.Coverage_pct == POLICY_COV)].iloc[0]
        cells.append(fyear(e['First_year_CIR_below4']))
        return " & ".join(cells) + " \\\\"
    body.append(unit_row(INDIA, bold=True))
    for region in REGION_ORDER:
        members = sorted([u for u in units if reg[u] == region], key=lambda u: -row(A, u, 0, 2026).CIR_mean)
        body.append("\\addlinespace[2pt]\\multicolumn{13}{l}{\\textit{" + region + "}} \\\\")
        body += [unit_row(u) for u in members]
    tab = tabular("l" + "c" * 12, hdr, body)
    caption = (f"State-level projections of cervical cancer incidence and mortality under a single-round HPV DNA testing programme at "
               f"{POLICY_COV}\\% coverage, by region, at three snapshot years. Observed rates are the 2012--19 means of the population-based cancer registry estimates used for "
               "calibration (ICMR-NINE); projected rates are per 100{,}000 women. Reductions are relative to no screening in the same year. "
               "First year CIR $<$ 4 is the first calendar year in which the mean CIR falls below 4 per 100{,}000 (the proxy for the WHO elimination threshold used in this study) at "
               f"{POLICY_COV}\\% coverage ({NEVER} = not reached by 2070). States are ordered within region by projected no-screening CIR in 2026. "
               "India row: dedicated national model; state rows: independently calibrated state models (1000 simulations each). The 2050 values also appear in Table~1 of the main text, "
               "beside those for VIA and Pap cytology. CIR=crude incidence rate. CMR=crude mortality rate.")
    note = ("Generated by hpvcc (paper tables). State models use the 18--85-year female population as denominator while the registries use all ages, "
            "so only rates are shown here; state-level case counts are not comparable with the India row.")
    return table_env(tab, caption, "tab:app_state_policy_hpv", landscape=True, size='\\scriptsize', note=note)

# ============================================================================= screening tests at policy coverage (Table 1, Table S11)
TEST_NOTE = ("Screening starts in 2026 and reaches 70\\% of women aged 30--65 years. The tests differ in the model only in their sensitivity "
             "for high-risk HPV infection ($H^+$): Pap cytology 0.70, visual inspection with acetic acid (VIA) 0.75, liquid-based cytology (LBC) 0.85 "
             "and HPV DNA testing 0.95. Rates are per 100{,}000 women; values are means of 1000 stochastic simulations with 95\\% intervals. ")
FIRST_NOTE = ("First year CIR $<$ 4 is the first calendar year in which the mean CIR falls below 4 per 100{,}000 (the proxy for the WHO "
              f"elimination threshold used in this study; {NEVER} = not reached by 2070). ")
MODEL_NOTE = "India row: dedicated national model; state rows: independently calibrated state models. CIR=crude incidence rate. CMR=crude mortality rate."


def _test_lookup(X):
    return {(r.Unit, r.Test_type): r for r in X.itertuples(index=False)}


def _test_ci(r, year):
    m, hw = getattr(r, f"CIR_{year}"), getattr(r, f"CIR_{year}_hw95")
    return fci(m, max(m - hw, 0.0), m + hw)


def _by_region(A, ncols, unit_row):
    """India first, then the states by region, ordered by no-screening CIR in 2026."""
    reg = region_of(A)
    body = [unit_row(INDIA, bold=True)]
    for region in REGION_ORDER:
        members = sorted([u for u in state_units(A) if reg[u] == region], key=lambda u: -row(A, u, 0, 2026).CIR_mean)
        body.append(f"\\addlinespace[2pt]\\multicolumn{{{ncols}}}{{l}}{{\\textit{{" + region + "}} \\\\")
        body += [unit_row(u) for u in members]
    return body


def t2_states(A, K, X):
    """Table 1: HPV DNA testing, VIA and Pap cytology at policy coverage, 2050."""
    tests = ["HPV DNA test", "VIA", "Pap cytology"]
    L = _test_lookup(X)
    def unit_row(u, bold=False):
        cells = [("\\textbf{" + tex_escape(u) + "}") if bold else tex_escape(u)]
        ko = K[K.Unit == u][['Observed_CIR', 'Observed_CMR']].mean()
        cells += [fnum(ko.Observed_CIR), fnum(ko.Observed_CMR), fnum(row(A, u, 0, 2050).CIR_mean)]
        cells += [_test_ci(L[u, t], 2050) for t in tests]
        cells += [fpct(L[u, t].CIR_reduction_pct_2050) for t in tests]
        cells += [fnum(L[u, t].CMR_2050) for t in tests]
        cells += [fyear(L[u, t].First_year_CIR_below4) for t in tests]
        return " & ".join(cells) + " \\\\"
    hdr = [f"State / Union Territory & \\multicolumn{{2}}{{c}}{{Observed 2012--19}} & CIR 2050, & \\multicolumn{{3}}{{c}}{{CIR 2050, {POLICY_COV}\\% coverage (95\\% CI)}} & "
           f"\\multicolumn{{3}}{{c}}{{CIR reduction 2050 (\\%)}} & \\multicolumn{{3}}{{c}}{{CMR 2050, {POLICY_COV}\\% coverage}} & \\multicolumn{{3}}{{c}}{{First year CIR $<$ 4}} \\\\",
           "\\cmidrule(lr){2-3}\\cmidrule(lr){5-7}\\cmidrule(lr){8-10}\\cmidrule(lr){11-13}\\cmidrule(lr){14-16}",
           " & CIR & CMR & no screening" + " & HPV DNA & VIA & Pap" * 4 + " \\\\"]
    tab = tabular("l" + "c" * 15, hdr, _by_region(A, 16, unit_row))
    caption = (f"State-level projections of cervical cancer incidence and mortality in 2050 under a single-round screening programme at {POLICY_COV}\\% coverage "
               "with HPV DNA testing, VIA or Pap cytology, by region. Observed rates are the 2012--19 means of the population-based cancer registry estimates used "
               "for calibration (ICMR-NINE). " + TEST_NOTE + "Reductions are relative to no screening in the same year. " + FIRST_NOTE +
               "States are ordered within region by projected no-screening CIR in 2026. HPV DNA testing at 2040 and 2060, and liquid-based cytology, are reported "
               "in the Supplementary Material (Appendix~K, Tables S10 and S11). " + MODEL_NOTE)
    return table_env(tab, caption, "tab:state_policy", landscape=True, size='\\scriptsize')


def s6_states_tests_other(A, X):
    """Table S11: what Table 1 and Table S10 do not show (VIA, Pap, LBC in 2040 and 2060; LBC in 2050)."""
    tests = ["VIA", "Pap cytology", "LBC"]
    L = _test_lookup(X)
    def unit_row(u, bold=False):
        cells = [("\\textbf{" + tex_escape(u) + "}") if bold else tex_escape(u)]
        cells += [fnum(L[u, t].CIR_2040) for t in tests] + [fnum(L[u, t].CIR_2060) for t in tests]
        lbc = L[u, "LBC"]
        cells += [_test_ci(lbc, 2050), fpct(lbc.CIR_reduction_pct_2050), fnum(lbc.CMR_2050), fyear(lbc.First_year_CIR_below4)]
        return " & ".join(cells) + " \\\\"
    hdr = [f"State / Union Territory & \\multicolumn{{3}}{{c}}{{CIR 2040, {POLICY_COV}\\% coverage}} & \\multicolumn{{3}}{{c}}{{CIR 2060, {POLICY_COV}\\% coverage}} & "
           f"\\multicolumn{{4}}{{c}}{{Liquid-based cytology, {POLICY_COV}\\% coverage}} \\\\",
           "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\\cmidrule(lr){8-11}",
           " & VIA & Pap & LBC & VIA & Pap & LBC & CIR 2050 (95\\% CI) & CIR reduction 2050 (\\%) & CMR 2050 & First year CIR $<$ 4 \\\\"]
    tab = tabular("l" + "c" * 10, hdr, _by_region(A, 11, unit_row))
    caption = (f"Screening tests other than HPV DNA testing at {POLICY_COV}\\% coverage, by state and Union Territory: projected crude incidence rate in 2040 and 2060 "
               "with VIA, Pap cytology and liquid-based cytology (LBC), and the 2050 results for LBC. This table complements Table~1 of the main text "
               "(HPV DNA testing, VIA and Pap cytology in 2050) and Table S10 (HPV DNA testing in 2040, 2050 and 2060). Screening starts in 2026; the tests differ "
               "in the model only in their sensitivity for high-risk HPV infection (Pap cytology 0.70, VIA 0.75, LBC 0.85). Rates are per 100{,}000 women "
               "(means of 1000 stochastic simulations, with the 95\\% interval where shown); the reduction is relative to no screening in the same year. "
               + FIRST_NOTE + MODEL_NOTE)
    return table_env(tab, caption, "tab:app_state_tests_other", landscape=True, size='\\scriptsize')


# ============================================================================= T3 elimination timeline
def t3_elimination(A, E, R):
    reg = region_of(A)
    units = state_units(A)
    covs = [25, 50, 70, 90, 100]
    hdr = ["State / Union Territory & \\multicolumn{5}{c}{First year with CIR $<$ 4 per 100{,}000, by coverage (\\%)} & "
           "\\multicolumn{3}{c}{Minimum coverage (\\%) for CIR $<$ 4 by} \\\\",
           "\\cmidrule(lr){2-6}\\cmidrule(lr){7-9}",
           " & 25 & 50 & 70 & 90 & 100 & 2040 & 2050 & 2070 \\\\"]
    def req(u, y):
        r = R[(R.Unit == u) & (R.Year == y)]
        if len(r) == 0:
            return NEVER
        r = r.iloc[0]
        return f"{int(r.Coverage_required_pct)}" if bool(r.Achievable_by_100pct) else "$>$100"
    def unit_row(u, bold=False):
        cells = [("\\textbf{" + tex_escape(u) + "}") if bold else tex_escape(u)]
        for c in covs:
            e = E[(E.Unit == u) & (E.Coverage_pct == c)].iloc[0]
            cells.append(fyear(e['First_year_CIR_below4']))
        cells += [req(u, y) for y in (2040, 2050, 2070)]
        return " & ".join(cells) + " \\\\"
    body = [unit_row(INDIA, bold=True)]
    for region in REGION_ORDER:
        members = sorted([u for u in units if reg[u] == region], key=lambda u: -row(A, u, 0, 2026).CIR_mean)
        body.append("\\addlinespace[2pt]\\multicolumn{9}{l}{\\textit{" + region + "}} \\\\")
        body += [unit_row(u) for u in members]
    tab = tabular("l" + "c" * 8, hdr, body)
    caption = ("Projected timing of cervical cancer elimination, assessed against a mean crude incidence rate below 4 per 100{,}000 women (a proxy for the WHO threshold, which is defined on age-standardised incidence), by state and Union Territory "
               f"under a single-round HPV DNA testing programme starting in 2026. Left: first calendar year below the threshold at each coverage ({NEVER} = not reached by 2070). "
               "Right: minimum coverage, in 5\\% steps, at which the threshold is met in the given year ($>$100 = not achievable at any coverage). "
               "Units whose no-screening CIR is already below 4 per 100{,}000 show 2026. India: dedicated national model; states: independently calibrated models.")
    note = "Generated by hpvcc (paper tables) from Elimination\\_Year and Coverage\\_Required."
    return table_env(tab, caption, "tab:elimination", landscape=False, size='\\footnotesize', note=note)

# ============================================================================= T4 regional summary (appendix)
def t4_regions(A):
    reg = region_of(A)
    units = state_units(A)
    hdr = ["Region & States & \\multicolumn{3}{c}{CIR, no screening} & \\multicolumn{3}{c}{CIR, 70\\% coverage} & "
           "\\multicolumn{3}{c}{CMR, no screening} & \\multicolumn{3}{c}{CMR, 70\\% coverage} & \\multicolumn{3}{c}{Units with CIR $<$ 4 at 70\\%} \\\\",
           "\\cmidrule(lr){3-5}\\cmidrule(lr){6-8}\\cmidrule(lr){9-11}\\cmidrule(lr){12-14}\\cmidrule(lr){15-17}",
           f" & (n) & {yh(YEARS)} & {yh(YEARS)} & {yh(YEARS)} & {yh(YEARS)} & {yh(YEARS)} \\\\"]
    def agg(members, cov, y, col):
        sub = A[(A.Unit.isin(members)) & (A.Coverage_pct == cov) & (A.Year == y)]
        return np.average(sub[col], weights=sub.Women_N)
    def nbelow(members, cov, y):
        sub = A[(A.Unit.isin(members)) & (A.Coverage_pct == cov) & (A.Year == y)]
        return int((sub.CIR_mean < 4).sum())
    body = []
    groups = [(r, [u for u in units if reg[u] == r]) for r in REGION_ORDER] + [("All 36 states/UTs (population-weighted)", units)]
    for name, members in groups:
        cells = [name, str(len(members))]
        cells += [fnum(agg(members, 0, y, 'CIR_mean')) for y in YEARS]
        cells += [fnum(agg(members, POLICY_COV, y, 'CIR_mean')) for y in YEARS]
        cells += [fnum(agg(members, 0, y, 'CMR_mean')) for y in YEARS]
        cells += [fnum(agg(members, POLICY_COV, y, 'CMR_mean')) for y in YEARS]
        cells += [f"{nbelow(members, POLICY_COV, y)}/{len(members)}" for y in YEARS]
        body.append(" & ".join(cells) + " \\\\")
    body.append("\\addlinespace[2pt]" + " & ".join(["India (dedicated national model)", "--"] +
                [fnum(row(A, INDIA, 0, y).CIR_mean) for y in YEARS] + [fnum(row(A, INDIA, POLICY_COV, y).CIR_mean) for y in YEARS] +
                [fnum(row(A, INDIA, 0, y).CMR_mean) for y in YEARS] + [fnum(row(A, INDIA, POLICY_COV, y).CMR_mean) for y in YEARS] +
                ["--"] * 3) + " \\\\")
    tab = tabular("l" + "c" * 16, hdr, body)
    caption = ("Regional summary of projected crude incidence and mortality rates (per 100{,}000 women) without screening and under a single-round HPV DNA testing "
               f"programme at 70\\% coverage, {YEARS[0]}--{YEARS[-1]}. Regional rates are population-weighted means of the independently calibrated state models (weights: modelled female "
               "population in each year). The last column counts the states and Union Territories in the region whose mean CIR is below 4 per 100{,}000 (proxy for the WHO elimination threshold). "
               "The final row is the dedicated national model used in the main text, which differs from the population-weighted sum of the state models.")
    note = "Generated by hpvcc (paper tables) from Annual\\_Long."
    return table_env(tab, caption, "tab:app_regions", landscape=True, size='\\scriptsize', note=note)

# ============================================================================= T5 incremental benefit (appendix)
def t5_incremental(A, C):
    steps = [(0, 25), (25, 50), (50, 70), (70, 90), (90, 100)]
    hdr = ["Coverage step & \\multicolumn{2}{c}{CIR reduction (per 100{,}000)} & \\multicolumn{2}{c}{CMR reduction (per 100{,}000)} & "
           "\\multicolumn{2}{c}{Additional cases averted 2026--2070} & \\multicolumn{2}{c}{Additional deaths averted 2026--2070} \\\\",
           "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}\\cmidrule(lr){8-9}",
           f"(\\%) & {yh(LATER)} & {yh(LATER)} & total & per 5 points of coverage & total & per 5 points of coverage \\\\"]
    body = []
    def cum(cov):
        c = C[(C.Unit == INDIA) & (C.Coverage_pct == cov) & (C['Horizon_year(2026..Y)'] == 2070)].iloc[0]
        return c.CumCases_averted_mean, c.CumDeaths_averted_mean
    for a, b in steps:
        cells = [f"{a} $\\rightarrow$ {b}"]
        for y in LATER:
            cells.append(fnum(row(A, INDIA, a, y).CIR_mean - row(A, INDIA, b, y).CIR_mean, 2))
        for y in LATER:
            cells.append(fnum(row(A, INDIA, a, y).CMR_mean - row(A, INDIA, b, y).CMR_mean, 2))
        ca, da = cum(a); cb, db = cum(b)
        n5 = (b - a) / 5
        cells += [fint(cb - ca), fint((cb - ca) / n5), fint(db - da), fint((db - da) / n5)]
        body.append(" & ".join(cells) + " \\\\")
    tab = tabular("l" + "c" * 8, hdr, body)
    caption = ("Incremental benefit of each step in programme coverage, India (dedicated national model). Rate reductions are the difference in projected crude rates "
               "between the lower and the higher coverage in the same year; cases and deaths averted are cumulative 2026--2070 and are also expressed per five percentage "
               "points of coverage to show the diminishing return above 70\\%.")
    note = "Generated by hpvcc (paper tables) from Annual\\_Long and Cumulative\\_Summary."
    return table_env(tab, caption, "tab:app_incremental", landscape=False, size='\\footnotesize', note=note)

# ============================================================================= T6 modality comparison (appendix)
def t6_modality(T):
    tt = T[T.Scenario.str.startswith('Test-type')]
    tests = [('Pap Cytology', 0.70), ('VIA', 0.75), ('LBC', 0.85), ('HPV DNA test', 0.95)]
    # match on H_pos to be robust to label spelling
    hdr = ["Screening test & $H^+$ & \\multicolumn{3}{c}{CIR (95\\% CI)} & \\multicolumn{3}{c}{CMR (95\\% CI)} & \\multicolumn{2}{c}{Averted 2026--2070 (vs no screening)} \\\\",
           "\\cmidrule(lr){3-5}\\cmidrule(lr){6-8}\\cmidrule(lr){9-10}",
           f" & & {yh(MODALITY_YEARS)} & {yh(MODALITY_YEARS)} & cases & deaths \\\\"]
    body = []
    cov = int(tt.Coverage_pct.iloc[0])
    for name, h in tests:
        sub = tt[np.isclose(tt.H_pos, h)]
        if len(sub) == 0:
            continue
        label = sub.Test_type.iloc[0]
        cells = [tex_escape(label), f"{h:.2f}"]
        for y in MODALITY_YEARS:
            r = sub[sub.Year == y].iloc[0]
            cells.append(fci(r.CIR_mean, r.CIR_lo95, r.CIR_hi95, 2))
        for y in MODALITY_YEARS:
            r = sub[sub.Year == y].iloc[0]
            cells.append(fci(r.CMR_mean, r.CMR_lo95, r.CMR_hi95, 2))
        s = sub[(sub.Year >= 2026) & (sub.Year <= 2070)]
        cells += [fint(s.Cases_averted_vs_noscreening.sum()), fint(s.Deaths_averted_vs_noscreening.sum())]
        body.append(" & ".join(cells) + " \\\\")
    tab = tabular("l" + "c" * 9, hdr, body)
    caption = (f"Projected national crude incidence and mortality rates (per 100{{,}}000 women) under four screening tests at a fixed coverage of {cov}\\%, "
               "India (dedicated national model). $H^+$ is the test sensitivity for high-risk HPV infection. Cases and deaths averted are summed over 2026--2070 from the annual "
               "ensemble means.")
    note = ("Generated by hpvcc (paper tables) from India\\_TestType (experiment india\\_test\\_type, the ensemble of the test-type figure). "
            "Intervals are mean $\\pm$ half the spread of the 2.5/97.5-percentile trajectories, as in the other tables.")
    return table_env(tab, caption, "tab:app_modality", landscape=False, size='\\footnotesize', note=note)

# ============================================================================= T7 calibration fit (appendix)
def t7_calibration(A, K, P):
    reg = region_of(A)
    units = state_units(A)
    g = K.groupby('Unit')
    hdr = ["State / Union Territory & \\multicolumn{3}{c}{CIR 2012--19 (mean)} & \\multicolumn{3}{c}{CMR 2012--19 (mean)} & "
           "\\multicolumn{4}{c}{Calibrated parameters} \\\\",
           "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\\cmidrule(lr){8-11}",
           " & Observed & Model & Diff.\\ (\\%) & Observed & Model & Diff.\\ (\\%) & $\\beta$ & $\\gamma$ ($\\times 10^{-3}$) & $r_I$ & $p$ \\\\"]
    Pi = P.set_index('Unit')
    def unit_row(u, bold=False):
        k = g.get_group(u)
        oc, mc, om, mm = k.Observed_CIR.mean(), k.Model_CIR.mean(), k.Observed_CMR.mean(), k.Model_CMR.mean()
        cells = [("\\textbf{" + tex_escape(u) + "}") if bold else tex_escape(u),
                 fnum(oc, 2), fnum(mc, 2), fnum(100 * (mc - oc) / oc, 1),
                 fnum(om, 2), fnum(mm, 2), fnum(100 * (mm - om) / om, 1)]
        p = Pi.loc[u]
        beta = p.beta_const if not pd.isna(p.beta_const) else p.get('beta_30_65', np.nan)
        cells += [fnum(beta, 3), fnum(1000 * p.gamma, 2), fnum(p.r_I, 3), fnum(p.p, 2)]
        return " & ".join(cells) + " \\\\"
    body = [unit_row(INDIA, bold=True)]
    for region in REGION_ORDER:
        members = sorted([u for u in units if reg[u] == region])
        body.append("\\addlinespace[2pt]\\multicolumn{11}{l}{\\textit{" + region + "}} \\\\")
        body += [unit_row(u) for u in members]
    tab = tabular("l" + "c" * 10, hdr, body)
    caption = ("Calibration fit and calibrated parameters for the national model and the 36 state and Union Territory models. Observed values are the registry-based "
               "estimates (ICMR-NINE) used for calibration; model values are the fitted no-screening trajectories, both averaged over the 2012--19 calibration window and "
               "expressed per 100{,}000 women. Diff.\\ is the relative difference of the model from the observed mean. $\\beta$ = HPV transmission rate; $\\gamma$ = progression "
               "rate from undetected infection to cancer; $r_I$ = clearance rate of undetected infection; $p$ = initial HPV prevalence used to initialise the infected compartment.")
    note = "Generated by hpvcc (paper tables) from Calibration\\_2012\\_2019 and Parameters."
    return table_env(tab, caption, "tab:app_calibration", landscape=False, size='\\scriptsize', note=note)


# ============================================================================= driver
TABLES = {
    "T1_national_coverage.tex": lambda f: t1_national(f["A"], f["C"], f["E"]),
    "T2_state_policy_coverage.tex": lambda f: t2_states(f["A"], f["K"], f["X"]),
    "T3_elimination_timeline.tex": lambda f: t3_elimination(f["A"], f["E"], f["R"]),
    "S1_regional_summary.tex": lambda f: t4_regions(f["A"]),
    "S2_incremental_benefit.tex": lambda f: t5_incremental(f["A"], f["C"]),
    "S3_modality_comparison.tex": lambda f: t6_modality(f["T"]),
    "S4_calibration_fit.tex": lambda f: t7_calibration(f["A"], f["K"], f["P"]),
    "S5_state_policy_hpv_dna.tex": lambda f: s5_states_hpv_dna(f["A"], f["E"], f["R"], f["K"]),
    "S6_state_tests_other.tex": lambda f: s6_states_tests_other(f["A"], f["X"]),
}
HEADER = "% Generated by hpvcc (hpvcc paper tables) -- do not edit by hand.\n"


def frames_from_csv(csvdir: Path) -> dict[str, pd.DataFrame]:
    """Read the master CSVs written by :meth:`hpvcc.paper.tables.MasterTables.save`."""
    csvdir = Path(csvdir)
    names = {"A": "Annual_Long", "E": "Elimination_Year", "R": "Coverage_Required", "C": "Cumulative_Summary",
             "K": "Calibration_2012_2019", "P": "Parameters", "T": "India_TestType", "X": "TestType_Elimination"}
    return {k: pd.read_csv(csvdir / f"{v}.csv") for k, v in names.items()}


def write_latex_tables(frames: dict[str, pd.DataFrame], out: Path | None = None, clean: bool = True) -> dict[str, Path]:
    """Write every table; ``frames`` maps A/E/R/C/K/P/T/X to the master tables."""
    global CLEAN
    CLEAN = clean
    out = ensure(Path(out) if out else TABLES_DIR / "latex")
    paths = {}
    for name, make in TABLES.items():
        p = out / name
        p.write_text(HEADER + make(frames), encoding="utf-8")
        paths[name] = p
    return paths
