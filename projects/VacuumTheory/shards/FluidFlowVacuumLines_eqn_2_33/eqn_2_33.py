from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_2_33__C_paralell_cap import eqn_2_33__C_paralell
from .eqn_2_33__arithmetic_sum_C_cap import eqn_2_33__arithmetic_sum_C

class FluidFlowVacuumLines:
    eqn_2_33__C_paralell = staticmethod(eqn_2_33__C_paralell)
    eqn_2_33__arithmetic_sum_C = staticmethod(eqn_2_33__arithmetic_sum_C)

    @kwasak_static
    def eqn_2_33(C_paralell=None, arithmetic_sum_C=None, **kwargs):
        return
