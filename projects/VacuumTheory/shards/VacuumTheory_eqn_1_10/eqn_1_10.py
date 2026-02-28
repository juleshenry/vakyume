from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_1_10__P_1_cap import eqn_1_10__P_1
from .eqn_1_10__P_2_cap import eqn_1_10__P_2
from .eqn_1_10__T_1_cap import eqn_1_10__T_1
from .eqn_1_10__T_2_cap import eqn_1_10__T_2
from .eqn_1_10__V_1_cap import eqn_1_10__V_1
from .eqn_1_10__V_2_cap import eqn_1_10__V_2

class VacuumTheory:
    eqn_1_10__P_1 = staticmethod(eqn_1_10__P_1)
    eqn_1_10__P_2 = staticmethod(eqn_1_10__P_2)
    eqn_1_10__T_1 = staticmethod(eqn_1_10__T_1)
    eqn_1_10__T_2 = staticmethod(eqn_1_10__T_2)
    eqn_1_10__V_1 = staticmethod(eqn_1_10__V_1)
    eqn_1_10__V_2 = staticmethod(eqn_1_10__V_2)

    @kwasak_static
    def eqn_1_10(P_1=None, P_2=None, T_1=None, T_2=None, V_1=None, V_2=None, **kwargs):
        return
