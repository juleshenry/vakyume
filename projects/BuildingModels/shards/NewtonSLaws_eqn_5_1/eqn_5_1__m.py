from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_5_1__m(self, F: float, a: float, **kwargs):
    # [.pyeqn] F = m * a
    result = []
    m = F / a
    result.append(m)
    return result
