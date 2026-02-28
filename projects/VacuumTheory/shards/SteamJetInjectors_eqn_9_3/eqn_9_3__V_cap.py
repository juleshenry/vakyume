from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_9_3__V(P_s: float, t_e: float, w_j: float, **kwargs):
    # [.pyeqn] t_e = (2.3 - 0.003 * P_s) * V / w_j
    result = []
    V = -1000.0*t_e*w_j/(3.0*P_s - 2300.0)
    result.append(V)
    return result
