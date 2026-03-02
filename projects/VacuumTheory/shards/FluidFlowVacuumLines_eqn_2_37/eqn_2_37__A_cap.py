from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_37__A(self, C: float, F_t: float, M: float, T: float, **kwargs):
    # [.pyeqn] C = 38.3 * (T * A * F_t / M) ** 0.5
    result = []
    A = 0.000681714375311032*C**2*M/(F_t*T)
    result.append(A)
    return result
