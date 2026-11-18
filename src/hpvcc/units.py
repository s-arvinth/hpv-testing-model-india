"""
Registry of modelled units: India (national model) and the 36 state/UT models.

Dadra and Nagar Haveli and Daman and Diu are modelled as one combined Union Territory
(merged in 2020); its inputs are the sums of the two former UTs (see
:mod:`hpvcc.data.states`).
"""
from __future__ import annotations

INDIA = "India"
COMBINED_UT = "Dadra and Nagar Haveli and Daman and Diu"

#: Six geographic regions and their member units (order is the paper-wide order).
REGIONS: dict[str, list[str]] = {
    "South": ["Karnataka", "Kerala", "Telangana", "Andhra Pradesh",
              "Tamil Nadu", "Puducherry", "Lakshadweep"],
    "North": ["Jammu and Kashmir", "Haryana", "Himachal Pradesh", "Uttarakhand",
              "Punjab", "Delhi", "Chandigarh", "Ladakh"],
    "Northeast": ["Assam", "Arunachal Pradesh", "Manipur", "Meghalaya",
                  "Mizoram", "Nagaland", "Sikkim", "Tripura"],
    "West": ["Gujarat", "Goa", "Maharashtra", "Rajasthan", COMBINED_UT],
    "East": ["Odisha", "Bihar", "West Bengal", "Jharkhand", "Andaman and Nicobar"],
    "Central": ["Madhya Pradesh", "Uttar Pradesh", "Chhattisgarh"],
}

#: All 36 state/UT model units, in region order.
STATES: list[str] = [s for members in REGIONS.values() for s in members]

#: Region of every state/UT model unit.
REGION_OF: dict[str, str] = {s: r for r, members in REGIONS.items() for s in members}

#: The two former UTs that make up the combined UT.
COMBINED_UT_PARTS = ("Daman and Diu", "Dadra and Nagar Haveli")


def slug(unit: str) -> str:
    """File-system name of a unit (``"Tamil Nadu" -> "Tamil_Nadu"``)."""
    return unit.replace(" ", "_")


def is_national(unit: str) -> bool:
    return unit == INDIA


def all_units() -> list[str]:
    """India followed by the 36 state/UT models."""
    return [INDIA, *STATES]


def resolve_unit(name: str) -> str:
    """Canonical unit name from a loose spelling (case, underscores, ``&``)."""
    key = name.strip().lower().replace("_", " ").replace("&", "and")
    for unit in all_units():
        if unit.lower() == key:
            return unit
    raise KeyError(f"unknown unit {name!r}; run `hpvcc units` for the list")


# --------------------------------------------------------------------------------------
# Name used for each unit in each raw data source (only entries that differ are listed).
# --------------------------------------------------------------------------------------

#: The NDMC population workbooks spell some names differently.
NDMC_TO_MODEL = {
    "Chattisgarh": "Chhattisgarh", "Orissa": "Odisha",
    "Jammu & Kashmir": "Jammu and Kashmir",
    "Andaman & Nicobar": "Andaman and Nicobar",
    "Dadra & Nagar Haveli": "Dadra and Nagar Haveli",
    "Daman & Diu": "Daman and Diu", "Puduchery": "Puducherry",
}
MODEL_TO_NDMC = {v: k for k, v in NDMC_TO_MODEL.items()}

#: State names in the population-projection workbook (Custom_AgeBands_All_States.xlsx).
MODEL_TO_POPPROJ = {
    "Andaman and Nicobar": "Andaman & Nicobar",
    "Chhattisgarh": "Chattisgarh",
    "Dadra and Nagar Haveli": "Dadra & Nagar Haveli",
    "Daman and Diu": "Daman & Diu",
    "Jammu and Kashmir": "Jammu & Kashmir",
    "Odisha": "Orissa",
    "Puducherry": "Puduchery",
}

#: IHME GBD reports some small units only as groups.
MODEL_TO_IHME = {
    "Jammu and Kashmir": "Jammu & Kashmir and Ladakh",
    "Ladakh": "Jammu & Kashmir and Ladakh",
    "Andaman and Nicobar": "Other Union Territories",
    "Chandigarh": "Other Union Territories",
    "Dadra and Nagar Haveli": "Other Union Territories",
    "Daman and Diu": "Other Union Territories",
    "Lakshadweep": "Other Union Territories",
    "Puducherry": "Other Union Territories",
}

#: Fixed population shares used to split the grouped IHME locations.
IHME_GROUP_SHARES = {
    "Jammu & Kashmir and Ladakh": {"Jammu and Kashmir": 0.979, "Ladakh": 0.021},
    "Other Union Territories": {
        "Andaman and Nicobar": 0.115, "Chandigarh": 0.303,
        "Dadra and Nagar Haveli": 0.100, "Daman and Diu": 0.070,
        "Lakshadweep": 0.021, "Puducherry": 0.391,
    },
}
