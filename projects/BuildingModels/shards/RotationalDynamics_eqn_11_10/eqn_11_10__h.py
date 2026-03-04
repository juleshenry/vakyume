from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_11_10__h(self, ICM: float, Ih: float, M: float, **kwargs):
    # [.pyeqn] Ih = ICM + M * h ** 2
    result = []
    h = sqrt((-ICM + Ih) / M)
    result.append(h)
    h = -sqrt(-(ICM - Ih) / M)
    result.append(h)
    return result
