from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_10_5__P_cap_1 import eqn_10_5__P_1
from .eqn_10_5__P_cap_2 import eqn_10_5__P_2
from .eqn_10_5__S_cap_p import eqn_10_5__S_p
from .eqn_10_5__V_cap import eqn_10_5__V
from .eqn_10_5__t import eqn_10_5__t


class LiquidRing:
    eqn_10_5__P_1 = eqn_10_5__P_1
    eqn_10_5__P_2 = eqn_10_5__P_2
    eqn_10_5__S_p = eqn_10_5__S_p
    eqn_10_5__V = eqn_10_5__V
    eqn_10_5__t = eqn_10_5__t

    @kwasak
    def eqn_10_5(self, P_1=None, P_2=None, S_p=None, V=None, t=None):
        return
