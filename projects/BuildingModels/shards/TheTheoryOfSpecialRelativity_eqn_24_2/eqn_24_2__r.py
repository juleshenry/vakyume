from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_24_2__r(self, FE: float, l: float, **kwargs):
    # [.pyeqn] FE = 2 * l / (2 * r)
    result = []
    r = l / FE
    result.append(r)
    return result
