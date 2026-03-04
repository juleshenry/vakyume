from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_3_6__H_cap_1 import eqn_3_6__H_1
from .eqn_3_6__H_cap_2 import eqn_3_6__H_2
from .eqn_3_6__P_cap import eqn_3_6__P
from .eqn_3_6__V_cap import eqn_3_6__V
from .eqn_3_6__V_cap_P_cap import eqn_3_6__V_P


class PressMgmt:
    eqn_3_6__H_1 = eqn_3_6__H_1
    eqn_3_6__H_2 = eqn_3_6__H_2
    eqn_3_6__P = eqn_3_6__P
    eqn_3_6__V = eqn_3_6__V
    eqn_3_6__V_P = eqn_3_6__V_P

    @kwasak
    def eqn_3_6(self, H_1=None, H_2=None, P=None, V=None, V_P=None):
        return
