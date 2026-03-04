from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_3_5__P import eqn_3_5__P
from .eqn_3_5__P_P import eqn_3_5__P_P
from .eqn_3_5__V import eqn_3_5__V
from .eqn_3_5__V_P import eqn_3_5__V_P

class PressMgmt:
    eqn_3_5__P = eqn_3_5__P
    eqn_3_5__P_P = eqn_3_5__P_P
    eqn_3_5__V = eqn_3_5__V
    eqn_3_5__V_P = eqn_3_5__V_P

    @kwasak
    def eqn_3_5(self, P=None, P_P=None, V=None, V_P=None):
        return
