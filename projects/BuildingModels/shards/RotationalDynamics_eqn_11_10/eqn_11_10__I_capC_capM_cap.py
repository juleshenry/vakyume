from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_11_10__ICM(self, Ih: float, M: float, h: float, **kwargs):
    # [.pyeqn] Ih = ICM + M * h ** 2
    result = []
    ICM = Ih - M * h**2
    result.append(ICM)
    return result
