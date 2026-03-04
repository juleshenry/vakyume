from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_3_5__P(self, P_P: float, V: float, V_P: float, **kwargs):
    # [.pyeqn] P_P = P * (V / V_P)
    result = []
    P = P_P*V_P/V
    result.append(P)
    return result
