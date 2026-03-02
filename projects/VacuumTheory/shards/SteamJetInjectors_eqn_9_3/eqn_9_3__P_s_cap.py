from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_9_3__P_s(self, V: float, t_e: float, w_j: float, **kwargs):
    # [.pyeqn] t_e = (2.3 - 0.003 * P_s) * V / w_j
    result = []
    P_s = 33.3333333333333*(23.0*V - 10.0*t_e*w_j)/V
    result.append(P_s)
    return result
