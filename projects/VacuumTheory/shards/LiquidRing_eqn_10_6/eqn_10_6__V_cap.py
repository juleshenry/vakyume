from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_10_6__V(P_1: float, P_2: float, S_a: float, t: float, **kwargs):
    # [.pyeqn] S_a = V / t * log(P_1 / P_2)
    result = []
    V = S_a*t/log(P_1/P_2)
    result.append(V)
    return result
