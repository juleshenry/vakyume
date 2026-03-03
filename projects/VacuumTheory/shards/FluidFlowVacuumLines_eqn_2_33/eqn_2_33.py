from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_2_33__C_cap_paralell import eqn_2_33__C_paralell
from .eqn_2_33__arithmetic_sum_C_cap import eqn_2_33__arithmetic_sum_C


class FluidFlowVacuumLines:
    eqn_2_33__C_paralell = eqn_2_33__C_paralell
    eqn_2_33__arithmetic_sum_C = eqn_2_33__arithmetic_sum_C

    @kwasak
    def eqn_2_33(self, C_paralell=None, arithmetic_sum_C=None):
        return
