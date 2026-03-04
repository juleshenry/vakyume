from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_10_6__P_cap_1 import eqn_10_6__P_1
from .eqn_10_6__P_cap_2 import eqn_10_6__P_2
from .eqn_10_6__S_cap_a import eqn_10_6__S_a
from .eqn_10_6__V_cap import eqn_10_6__V
from .eqn_10_6__t import eqn_10_6__t

class LiquidRing:
    eqn_10_6__P_1 = eqn_10_6__P_1
    eqn_10_6__P_2 = eqn_10_6__P_2
    eqn_10_6__S_a = eqn_10_6__S_a
    eqn_10_6__V = eqn_10_6__V
    eqn_10_6__t = eqn_10_6__t

    @kwasak
    def eqn_10_6(self, P_1=None, P_2=None, S_a=None, V=None, t=None):
        return
