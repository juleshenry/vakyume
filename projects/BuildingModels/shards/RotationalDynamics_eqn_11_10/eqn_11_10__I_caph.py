from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_11_10__Ih(self, ICM: float, M: float, h: float, **kwargs):
    # [.pyeqn] Ih = ICM + M * h ** 2
    result = []
    Ih = ICM + M * h**2
    result.append(Ih)
    return result
