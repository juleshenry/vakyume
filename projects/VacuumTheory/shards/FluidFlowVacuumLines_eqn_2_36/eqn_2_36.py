from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_2_36__C_cap import eqn_2_36__C
from .eqn_2_36__C_0_cap import eqn_2_36__C_0
from .eqn_2_36__F_t_cap import eqn_2_36__F_t

class FluidFlowVacuumLines:
    eqn_2_36__C = eqn_2_36__C
    eqn_2_36__C_0 = eqn_2_36__C_0
    eqn_2_36__F_t = eqn_2_36__F_t

    @kwasak_static
    def eqn_2_36(self, C=None, C_0=None, F_t=None):
        return
