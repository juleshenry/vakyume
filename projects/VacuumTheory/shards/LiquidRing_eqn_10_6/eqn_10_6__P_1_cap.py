from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_10_6__P_1(P_2: float, S_a: float, V: float, t: float, **kwargs):
    # [.pyeqn] S_a = V / t * log(P_1 / P_2)
    result = []
    P_1 = P_2*exp(S_a*t/V)
    result.append(P_1)
    return result
