from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_3_5__P_cap import eqn_3_5__P
from .eqn_3_5__P_P_cap import eqn_3_5__P_P
from .eqn_3_5__V_cap import eqn_3_5__V
from .eqn_3_5__V_P_cap import eqn_3_5__V_P

class PressMgmt:
    eqn_3_5__P = staticmethod(eqn_3_5__P)
    eqn_3_5__P_P = staticmethod(eqn_3_5__P_P)
    eqn_3_5__V = staticmethod(eqn_3_5__V)
    eqn_3_5__V_P = staticmethod(eqn_3_5__V_P)

    @kwasak_static
    def eqn_3_5(P=None, P_P=None, V=None, V_P=None, **kwargs):
        return
