# Data

The model needs, for each unit, projected women by age band, a rate denominator, a
recruitment inflow and an initial state in 2002. These **model inputs** are in
[`processed/units/`](processed/units) (one JSON file per unit) and are built from the
**source files** in [`raw/`](raw) by

```bash
hpvcc prepare-inputs
```

The test suite checks that rebuilding reproduces the stored inputs exactly.

## Source files (`raw/`)

| File | Content | Source | Used for |
|---|---|---|---|
| `popproj_female_agebands_states.xlsx` | women by state and model age band (18–29, 30–65, 66+), 1991–2100 | female population projections, <https://github.com/bisect-group/-india-female-pop-projections> | state/UT populations, recruitment and denominators |
| `ndmc_india_female_agebands.xlsx` | women aged 18–29, 30–65, 66+ in India, 2012–2070 | NDMC population projections | national-model populations and recruitment |
| `ndmc_female_model_agebands_2012-2070.xlsx` | the same three bands by state (`States_Female_Long`) and India (`India_Female`) | NDMC population projections | state initial infections I(2002) |
| `ndmc_statewise_population.xlsx` | all-age women by state, 2012–2036 (`Female`) | NDMC population projections | national rate denominator |
| `ihme_gbd2023_cervical_cancer_prevalence.csv` | prevalent cervical cancers by state, 1990–2023 | IHME Global Burden of Disease 2023 | state initial cancers C(2002) |
| `state_initial_hpv_prevalence.csv` | HPV prevalence (%) per state/UT and each unit's share (%) of the national population | literature review (paper, Appendix Table S2) | initial infections I(2002) |

How each input is derived (back-casting to 2002, the ten-year recruitment lead, the
combined UT of Dadra and Nagar Haveli and Daman and Diu, IHME groups split by population
share) is documented in [`src/hpvcc/data/states.py`](../src/hpvcc/data/states.py) and
[`src/hpvcc/data/national.py`](../src/hpvcc/data/national.py), and summarised in
[`docs/model.md`](../docs/model.md#inputs).

## Observed rates (`registry/`)

[`registry/observed_rates_2012_2019.csv`](registry/observed_rates_2012_2019.csv) holds the
values the models were fitted to: crude cervical-cancer incidence (`CIR`) and mortality
(`CMR`) rates per 100,000 women, for India and each of the 36 state/UT units, 2012–2019,
from the population-based cancer registry (PBCR) estimates of ICMR's National Centre for
Disease Informatics and Research (NCDIR) / National Institute of Cancer Epidemiology and
Research (NINE). Tamil Nadu, Mizoram, Manipur, Sikkim and Tripura are directly observed
by a PBCR; the other units are PBCR-based estimates.

They are not model inputs: the models run from their fitted parameters (`params/`). The
rates are used to show the fit — `hpvcc fit UNIT`, [`hpvcc.observed`](../src/hpvcc/observed.py),
Figs 1, 2 and S5–S10, Table 1 (observed columns) and Table S9. Case and death counts and
age-group data are not included, and neither is the calibration code.
