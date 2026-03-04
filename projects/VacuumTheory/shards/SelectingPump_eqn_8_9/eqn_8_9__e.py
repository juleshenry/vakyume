from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_8_9__e(self, E_j: float, E_m: float, r: float, s: float, **kwargs):
    # [.pyeqn] r = 2.93 * (E_j * e) / (E_m * s)
    result = []
    e = 0.341296928327645 * E_m * r * s / E_j
    result.append(e)
    return result
