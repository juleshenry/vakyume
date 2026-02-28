from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_10_6__P_1_cap import eqn_10_6__P_1
from .eqn_10_6__P_2_cap import eqn_10_6__P_2
from .eqn_10_6__S_a_cap import eqn_10_6__S_a
from .eqn_10_6__V_cap import eqn_10_6__V
from .eqn_10_6__t import eqn_10_6__t

class LiquidRing:
    eqn_10_6__P_1 = staticmethod(eqn_10_6__P_1)
    eqn_10_6__P_2 = staticmethod(eqn_10_6__P_2)
    eqn_10_6__S_a = staticmethod(eqn_10_6__S_a)
    eqn_10_6__V = staticmethod(eqn_10_6__V)
    eqn_10_6__t = staticmethod(eqn_10_6__t)

    @kwasak_static
    def eqn_10_6(P_1=None, P_2=None, S_a=None, V=None, t=None, **kwargs):
        return
