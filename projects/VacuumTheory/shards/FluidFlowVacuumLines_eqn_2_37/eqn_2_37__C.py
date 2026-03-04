from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_2_37__C(self, A: float, F_t: float, M: float, T: float, **kwargs):
    # [.pyeqn] C = 38.3 * (T * A * F_t / M) ** 0.5
    result = []
    C = 38.3 * sqrt(A * F_t * T / M)
    result.append(C)
    return result
