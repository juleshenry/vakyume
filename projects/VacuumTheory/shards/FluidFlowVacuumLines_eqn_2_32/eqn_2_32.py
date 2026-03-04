from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_2_32__C_cap_series import eqn_2_32__C_series
from .eqn_2_32__geometric_sum_C_cap import eqn_2_32__geometric_sum_C

class FluidFlowVacuumLines:
    eqn_2_32__C_series = eqn_2_32__C_series
    eqn_2_32__geometric_sum_C = eqn_2_32__geometric_sum_C

    @kwasak
    def eqn_2_32(self, C_series=None, geometric_sum_C=None):
        return
