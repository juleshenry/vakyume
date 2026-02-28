from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_2_8__M_cap import eqn_2_8__M
from .eqn_2_8__P_c_cap import eqn_2_8__P_c
from .eqn_2_8__T_c_cap import eqn_2_8__T_c
from .eqn_2_8__mu_c import eqn_2_8__mu_c

class FluidFlowVacuumLines:
    eqn_2_8__M = staticmethod(eqn_2_8__M)
    eqn_2_8__P_c = staticmethod(eqn_2_8__P_c)
    eqn_2_8__T_c = staticmethod(eqn_2_8__T_c)
    eqn_2_8__mu_c = staticmethod(eqn_2_8__mu_c)

    @kwasak_static
    def eqn_2_8(M=None, P_c=None, T_c=None, mu_c=None, **kwargs):
        return
