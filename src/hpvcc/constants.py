"""Model-wide constants shared by every unit (national and state/UT models)."""

#: Compartment names, in the order used by the ODE state vector.
COMPARTMENTS = ("S18", "S30", "S61", "S_R", "H_T", "I", "C", "R_H", "C_S", "C_D")

#: Index of each compartment in the state vector.
IDX = {name: i for i, name in enumerate(COMPARTMENTS)}

#: Human-readable compartment names (used in figure titles).
COMPARTMENT_LABELS = {
    "S18": "Susceptible (18–29)",
    "S30": "Susceptible (30–65)",
    "S61": "Susceptible (66–85)",
    "S_R": "Reinfection-Susceptible",
    "H_T": "Known Infected",
    "I": "Unknown Infected",
    "C": "Cervical Cancer",
    "R_H": "HPV Recovered",
    "C_S": "Cancer Survivors",
    "C_D": "Cumulative Cancer Deaths",
}

#: Mathematical symbol of each compartment (LaTeX / matplotlib mathtext).
COMPARTMENT_TEX = {
    "S18": "S_{18}", "S30": "S_{30}", "S61": "S_{61}", "S_R": "S_R", "H_T": "H_T",
    "I": "I", "C": "C", "R_H": "R_H", "C_S": "C_S", "C_D": "C_D",
}

SIM_START = 2002        #: first simulated year
CALIB_START = 2012      #: first year of the calibration window (also the population back-cast window)
CALIB_END = 2019        #: last year of the calibration window
SIM_END = 2070          #: projection horizon
SCREENING_START = 2026  #: first year of the screening programme
N_BACK = CALIB_START - SIM_START   #: years simulated before the calibration window

MU_NATURAL = 0.01388    #: natural (background) mortality rate, all compartments, per year
WHO_THRESHOLD = 4.0     #: elimination proxy: crude incidence rate < 4 per 100,000 women
H_POS_HPV_DNA = 0.95    #: sensitivity of HPV-DNA testing (default screening test)

#: Screening tests and their sensitivity H+ (probability that an infected, screened woman tests positive).
SCREENING_TESTS = {
    "Pap cytology": 0.70,
    "VIA": 0.75,
    "LBC": 0.85,
    "HPV DNA test": 0.95,
}

#: Scale factor for rates: model populations are in thousands, rates are per 100,000.
PER_100K = 1e5

#: Short names accepted by the command line for each screening test.
TEST_ALIASES = {"pap": "Pap cytology", "via": "VIA", "lbc": "LBC", "hpv-dna": "HPV DNA test", "hpv": "HPV DNA test"}
