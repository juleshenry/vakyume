from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_37__F_t(self, A: float, C: float, M: float, T: float, **kwargs):
    # [.pyeqn] C = 38.3 * (T * A * F_t / M) ** 0.5
    result = []
    F_t = 0.000681714375311032*C**2*M/(A*T)
    result.append(F_t)
    return result
