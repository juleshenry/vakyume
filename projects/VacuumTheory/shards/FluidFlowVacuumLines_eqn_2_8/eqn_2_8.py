from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_2_8__M import eqn_2_8__M
from .eqn_2_8__P_c import eqn_2_8__P_c
from .eqn_2_8__T_c import eqn_2_8__T_c
from .eqn_2_8__mu_c import eqn_2_8__mu_c

class FluidFlowVacuumLines:
    eqn_2_8__M = eqn_2_8__M
    eqn_2_8__P_c = eqn_2_8__P_c
    eqn_2_8__T_c = eqn_2_8__T_c
    eqn_2_8__mu_c = eqn_2_8__mu_c

    @kwasak
    def eqn_2_8(self, M=None, P_c=None, T_c=None, mu_c=None):
        """
        M:= mol. weight
        T_c:= critical temp, K
        P_c:= critical pressure, atm
        """
        return
