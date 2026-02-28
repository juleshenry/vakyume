from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_2_32__C_series_cap import eqn_2_32__C_series
from .eqn_2_32__geometric_sum_C_cap import eqn_2_32__geometric_sum_C

class FluidFlowVacuumLines:
    eqn_2_32__C_series = staticmethod(eqn_2_32__C_series)
    eqn_2_32__geometric_sum_C = staticmethod(eqn_2_32__geometric_sum_C)

    @kwasak_static
    def eqn_2_32(C_series=None, geometric_sum_C=None, **kwargs):
        return
