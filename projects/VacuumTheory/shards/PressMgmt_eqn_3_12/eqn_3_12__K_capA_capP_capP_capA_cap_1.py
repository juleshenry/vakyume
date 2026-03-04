from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_3_12__KAPPA_1(self, H_2: float, P: float, **kwargs):
    # [.pyeqn] P = KAPPA_1 * H_2 ** 2
    result = []
    KAPPA_1 = P / H_2**2
    result.append(KAPPA_1)
    return result
