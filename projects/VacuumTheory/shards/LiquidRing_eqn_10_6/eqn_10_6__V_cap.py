from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_10_6__V(self, P_1: float, P_2: float, S_a: float, t: float, **kwargs):
    # [.pyeqn] S_a = V / t * log(P_1 / P_2)
    result = []
    V = S_a*t/log(P_1/P_2)
    result.append(V)
    return result
