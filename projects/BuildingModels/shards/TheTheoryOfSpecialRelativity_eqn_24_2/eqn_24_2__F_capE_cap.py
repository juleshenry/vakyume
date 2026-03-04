from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_24_2__FE(self, l: float, r: float, **kwargs):
    # [.pyeqn] FE = 2 * l / (2 * r)
    result = []
    FE = l / r
    result.append(FE)
    return result
