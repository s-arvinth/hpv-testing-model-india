"""
hpvcc — a compartmental model of HPV infection, HPV testing and cervical cancer in women,
for India (national model) and its 36 states and Union Territories.

Ten compartments (susceptible women in three age bands, re-susceptible, screen-detected
and undetected HPV infection, cervical cancer, recovered, cancer survivors and cumulative
cancer deaths), one fitted parameter set per unit, projections 2002–2070 with optional
screening (coverage α, test sensitivity H⁺) from 2026.

Quick start
-----------
>>> from hpvcc import project
>>> proj = project("India", coverage=[0, 0.7])
>>> proj.summary()

Main building blocks
--------------------
* :class:`~hpvcc.parameters.ModelParameters`, :func:`~hpvcc.parameters.load_parameters`
* :class:`~hpvcc.inputs.UnitInputs` — population, recruitment, initial conditions
* :class:`~hpvcc.model.HPVCervicalCancerModel` — the ODE system
* :mod:`hpvcc.simulation` — scenarios and stochastic ensembles
* :func:`~hpvcc.projection.project` — scenario → table of CIR / CMR
* :mod:`hpvcc.uncertainty` — parameter-uncertainty band of the national model
* :mod:`hpvcc.paper` — reproduction of the paper's simulations, tables and figures
"""

__version__ = "1.0.0"

from .model import HPVCervicalCancerModel  # noqa: E402
from .parameters import ModelParameters, load_parameters  # noqa: E402
from .projection import Projection, project  # noqa: E402
from .simulation import Scenario  # noqa: E402

__all__ = ["HPVCervicalCancerModel", "ModelParameters", "load_parameters", "Projection", "project",
           "Scenario", "__version__"]
