from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_10_5__P_2(self, P_1: float, S_p: float, V: float, t: float, **kwargs):
    # [.pyeqn] t = V / S_p * log(P_1 / P_2)
    result = []
    P_2 = P_1*exp(-S_p*t/V)
    result.append(P_2)
    return result
