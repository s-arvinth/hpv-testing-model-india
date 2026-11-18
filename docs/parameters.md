# Parameters

Rates are per year; populations are in thousands of women. Field names are those of
[`ModelParameters`](../src/hpvcc/parameters.py); CSV column names in brackets where they
differ.

| Symbol | Field | Meaning | Fitted? |
|:---:|---|---|---|
| <i>p</i> | `p` | scaling of the recruitment inflow <i>Λ</i>(<i>t</i>) into <i>S</i><sub>18</sub> | per unit |
| <i>K</i><sub>1</sub> | `K1` | ageing 18–29 → 30–65 | per unit |
| <i>K</i><sub>2</sub> | `K2` | ageing 30–65 → 66–85 | per unit |
| <i>β</i> | `beta` (`beta_const`) | HPV transmission rate, the same in the three age bands | per unit |
| <i>β</i><sub>R</sub> | `beta_reinfection` (`beta_Reinfec`) | transmission rate for re-susceptible women | per unit |
| <i>γ</i> | `gamma` | progression of undetected infection to cancer | per unit |
| <i>θ</i> | — | progression of screen-detected infection to cancer | **sampled**: U(0, <i>γ</i>), redrawn every year |
| <i>r</i><sub>I</sub> | `r_I` | clearance of undetected infection | per unit |
| <i>r</i><sub>H</sub> | — | clearance / treatment of screen-detected infection | **sampled**: U(<i>r</i><sub>I</sub>, 1), redrawn every year |
| <i>r</i><sub>C</sub> | `r_C` | cancer → survivor | per unit |
| <i>μ</i><sub>C</sub> | `mu_C` | cancer-mortality coefficient (deaths = <i>μ</i><sub>C</sub> <i>C</i><sup><i>n</i><sub>PL</sub></sup>) | per unit |
| <i>n</i><sub>PL</sub> | `n_PL` | exponent of the cancer-mortality power law | per unit |
| <i>ω</i> | `omega` | relapse of survivors to cancer | per unit |
| <i>ε</i> | `epsilon` | loss of post-clearance immunity (<i>R</i><sub>H</sub> → <i>S</i><sub>R</sub>) | per unit |
| <i>μ</i> | `mu` | natural mortality | fixed, 0.01388 |
| <i>H</i><sup>+</sup> | `H_pos` | screening-test sensitivity | scenario: HPV-DNA 0.95, LBC 0.85, VIA 0.75, Pap 0.70 |

## Files

| File | Content |
|---|---|
| [`params/states_final.csv`](../params/states_final.csv) | 36 state/UT models (one row each, column `State`) |
| [`params/india_final.csv`](../params/india_final.csv) | national model (also lists `mu`, `H_pos`) |
| [`params/india_uncertainty_sets.csv`](../params/india_uncertainty_sets.csv) | 500 accepted national parameter sets (parameter uncertainty), best first |

Values are stored at full double precision; read them with
`pandas.read_csv(..., float_precision="round_trip")`, as the package does. State rows have
no `mu` or `H_pos` columns and use μ = 0.01388, H⁺ = 0.95. No file holds a value of θ or
r_H: these two rates of the screened pathway are sampled in every simulation
([model](model.md#stochastic-screened-pathway)).

## Using other values

```python
from hpvcc import load_parameters, project, HPVCervicalCancerModel

p = load_parameters("Bihar")                                   # fitted set
project("Bihar", coverage=0.7, changes={"beta": 0.5})          # change some values
project("Bihar", coverage=0.7, params=p.replace(gamma=0.001))  # a full set
project("Bihar", coverage=0.7, params_file="my_params.csv")    # your own CSV
model = HPVCervicalCancerModel.for_unit("Bihar").with_params(r_I=0.6)
```

```bash
hpvcc run Bihar --coverage 0.7 --set beta=0.5 gamma=0.001
hpvcc run Bihar --coverage 0.7 --params-file my_params.csv
```

A parameter CSV of your own may hold one row (used for any unit) or one row per unit with
a `State` (or `Unit`) column; columns follow `params/*.csv` or the field names above.

## Fitted values

Rounded to four significant digits (full precision in `params/`).

| Unit | Region | β | β_R | γ | r_I | r_C | μ_C | n_PL | ω | ε | p | K1 | K2 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| **India** | national | 0.5 | 0.5 | 0.00058 | 0.4785 | 0.1 | 22 | 0.1 | 0.01 | 0.07 | 0.082 | 0.06761 | 0.00836 |
| Andhra Pradesh | South | 0.82 | 0.8 | 0.0024 | 0.7631 | 0.2 | 0.09 | 1.45 | 0.1 | 1.267 | 0.22 | 0.17 | 0.014 |
| Karnataka | South | 0.53 | 0.5 | 0.00075 | 0.435 | 0.7915 | 0.9731 | 0.9667 | 0.24 | 0.5653 | 0.21 | 0.18 | 0.015 |
| Kerala | South | 0.283 | 0.1511 | 0.0009834 | 0.2311 | 0.9767 | 1.395 | 1.342 | 0.1805 | 0.4856 | 0.17 | 0.15 | 0.015 |
| Lakshadweep | South | 0.1809 | 0.3616 | 0.002633 | 0.1822 | 0.9951 | 14.12 | 1.377 | 0.09475 | 1.926 | 0.1363 | 0.115 | 0.01097 |
| Puducherry | South | 0.13 | 0.15 | 0.0017 | 0.105 | 0.951 | 1.261 | 1.034 | 0.1135 | 1.265 | 0.15 | 0.1236 | 0.01051 |
| Tamil Nadu | South | 0.2 | 0.2 | 0.0018 | 0.16 | 0.2892 | 1.253 | 0.4192 | 0 | 0.3321 | 0.17 | 0.15 | 0.01 |
| Telangana | South | 0.22 | 0.21 | 0.0026 | 0.1849 | 0.2508 | 0.3085 | 1.024 | 0.005504 | 0.4747 | 0.14 | 0.11 | 0.01 |
| Chandigarh | North | 0.23 | 0.23 | 0.0003 | 0.2 | 0.7969 | 4.533 | 1.361 | 0.0432 | 0.5646 | 0.1 | 0.07 | 0.007 |
| Delhi | North | 0.1004 | 0.8646 | 0.002139 | 0.1936 | 0.29 | 0.33 | 1.046 | 0.1772 | 1.55 | 0.09934 | 0.08252 | 0.008209 |
| Haryana | North | 0.64 | 0.6 | 0.00088 | 0.57 | 0.6291 | 25 | 5 | 0.1853 | 1.221 | 0.15 | 0.11 | 0.012 |
| Himachal Pradesh | North | 0.1373 | 0.1076 | 0.0009978 | 0.1119 | 0.7339 | 1.621 | 1.313 | 0.06447 | 0.785 | 0.11 | 0.09 | 0.007594 |
| Jammu and Kashmir | North | 0.101 | 0.55 | 0.00011 | 0.1404 | 0.932 | 4.454 | 1.376 | 0.005153 | 1.226 | 0.09 | 0.06832 | 0.008207 |
| Ladakh | North | 0.1271 | 0.1092 | 0.0001253 | 0.1005 | 0.5313 | 31.07 | 0.782 | 0.9712 | 0.2815 | 0.106 | 0.08306 | 0.009396 |
| Punjab | North | 0.235 | 0.1035 | 0.00045 | 0.1676 | 0.7265 | 0.9563 | 1.179 | 0.2974 | 0.6488 | 0.175 | 0.14 | 0.015 |
| Uttarakhand | North | 0.198 | 0.1173 | 0.0005542 | 0.1342 | 0.9843 | 2.473 | 1.396 | 0.9958 | 0.04088 | 0.13 | 0.1 | 0.01 |
| Arunachal Pradesh | Northeast | 0.12 | 0.12 | 0.0005181 | 0.1051 | 0.9219 | 6.894 | 1.395 | 0.09959 | 0.6559 | 0.07 | 0.05 | 0.003689 |
| Assam | Northeast | 0.226 | 0.102 | 0.001178 | 0.1722 | 0.9044 | 1.117 | 1.07 | 0.1571 | 0.7419 | 0.1 | 0.07 | 0.0064 |
| Manipur | Northeast | 0.135 | 0.13 | 0.0005193 | 0.1064 | 0.5006 | 0.8084 | 1.095 | 0.9632 | 1.866 | 0.09 | 0.05938 | 0.006666 |
| Meghalaya | Northeast | 0.16 | 0.135 | 0.00016 | 0.115 | 0.98 | 7.937 | 1.493 | 0.996 | 1.165 | 0.09 | 0.06 | 0.006 |
| Mizoram | Northeast | 0.14 | 0.15 | 0.001677 | 0.1176 | 0.8052 | 4.742 | 1.498 | 0.05031 | 1.835 | 0.085 | 0.065 | 0.0065 |
| Nagaland | Northeast | 0.15 | 0.15 | 0.001739 | 0.127 | 0.7611 | 1.318 | 1.133 | 0.1469 | 0.4594 | 0.09 | 0.07 | 0.006 |
| Sikkim | Northeast | 0.2 | 0.2 | 0.0005168 | 0.18 | 0.7927 | 10.5 | 1.5 | 0.937 | 1.852 | 0.1 | 0.07 | 0.01 |
| Tripura | Northeast | 0.135 | 0.135 | 0.0009466 | 0.112 | 0.302 | 0.9223 | 1.492 | 0.04474 | 0.8762 | 0.08 | 0.06 | 0.007 |
| Dadra and Nagar Haveli and Daman and Diu | West | 0.212 | 0.18 | 0.0015 | 0.1738 | 0.5525 | 2.5 | 1.303 | 0.4495 | 0.2906 | 0.1301 | 0.1169 | 0.01286 |
| Goa | West | 0.25 | 0.25 | 0.00095 | 0.25 | 0.5866 | 1.54 | 1.244 | 0.2931 | 1.819 | 0.12 | 0.1 | 0.009 |
| Gujarat | West | 0.5 | 0.45 | 0.0008 | 0.4268 | 0.1087 | 0.65 | 0.2333 | 0.1364 | 0.2 | 0.17 | 0.14 | 0.015 |
| Maharashtra | West | 0.1003 | 0.8794 | 0.001412 | 0.1845 | 0.1564 | 0.4059 | 0.6662 | 0.08442 | 1.379 | 0.1086 | 0.09101 | 0.009464 |
| Rajasthan | West | 0.65 | 0.63 | 0.0016 | 0.59 | 0.006138 | 0.25 | 0.48 | 0.2451 | 1.109 | 0.16 | 0.13 | 0.015 |
| Andaman and Nicobar | East | 0.129 | 0.1095 | 0.001955 | 0.1109 | 0.5221 | 4.299 | 1.49 | 0.1089 | 1.024 | 0.1245 | 0.1037 | 0.008608 |
| Bihar | East | 0.68 | 0.56 | 0.0008 | 0.56 | 0.8355 | 0.5 | 2.5 | 0.1549 | 1.1 | 0.16 | 0.12 | 0.016 |
| Jharkhand | East | 0.361 | 0.3 | 0.00092 | 0.291 | 0.3025 | 0.3 | 1.2 | 0.3077 | 0.262 | 0.13 | 0.1 | 0.009 |
| Odisha | East | 0.155 | 0.09 | 0.00045 | 0.1075 | 0.9297 | 1.109 | 1.043 | 0.1734 | 0.9576 | 0.16 | 0.13 | 0.015 |
| West Bengal | East | 0.82 | 0.71 | 0.00205 | 0.7236 | 0.9842 | 0.905 | 1.332 | 0.088 | 0.6966 | 0.19 | 0.15 | 0.015 |
| Chhattisgarh | Central | 0.32 | 0.28 | 0.00065 | 0.27 | 0.95 | 4 | 2.5 | 0.1569 | 0.9627 | 0.15 | 0.12 | 0.012 |
| Madhya Pradesh | Central | 0.33 | 0.22 | 0.0007 | 0.263 | 0.6573 | 0.7992 | 0.9818 | 0.12 | 1.996 | 0.14 | 0.12 | 0.0085 |
| Uttar Pradesh | Central | 0.47 | 0.37 | 0.00125 | 0.4 | 0.818 | 0.9329 | 1.034 | 0.1508 | 1.94 | 0.11 | 0.08 | 0.008 |
