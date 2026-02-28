from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_3_6__H_1_cap import eqn_3_6__H_1
from .eqn_3_6__H_2_cap import eqn_3_6__H_2
from .eqn_3_6__P_cap import eqn_3_6__P
from .eqn_3_6__V_cap import eqn_3_6__V
from .eqn_3_6__V_P_cap import eqn_3_6__V_P

class PressMgmt:
    eqn_3_6__H_1 = staticmethod(eqn_3_6__H_1)
    eqn_3_6__H_2 = staticmethod(eqn_3_6__H_2)
    eqn_3_6__P = staticmethod(eqn_3_6__P)
    eqn_3_6__V = staticmethod(eqn_3_6__V)
    eqn_3_6__V_P = staticmethod(eqn_3_6__V_P)

    @kwasak_static
    def eqn_3_6(H_1=None, H_2=None, P=None, V=None, V_P=None, **kwargs):
        return
