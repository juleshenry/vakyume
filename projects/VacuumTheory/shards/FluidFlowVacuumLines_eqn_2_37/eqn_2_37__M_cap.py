from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_37__M(self, A: float, C: float, F_t: float, T: float, **kwargs):
    # [.pyeqn] C = 38.3 * (T * A * F_t / M) ** 0.5
    result = []
    M = 1466.89*A*F_t*T/C**2
    result.append(M)
    return result
