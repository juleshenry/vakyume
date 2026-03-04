from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_8_9__r(self, E_j: float, E_m: float, e: float, s: float, **kwargs):
    # [.pyeqn] r = 2.93 * (E_j * e) / (E_m * s)
    result = []
    r = 2.93 * E_j * e / (E_m * s)
    result.append(r)
    return result
