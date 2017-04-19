"""
The :mod:`steapy` package implements a boundary element method for steady inviscid
fluid flows.
"""

from .airfoil_generate import *
from .panels import *
from .freestream import *
from .integral import *
from .source_contribution import *
from .vortex_contribution import *
from .kutta_condition import *
from .singularity_matrix import *
from .build_rhs import *
from .tangential_vel import *
from .pressure_coeff import *
from .velocity_field import *

__version__ = "0.1"
