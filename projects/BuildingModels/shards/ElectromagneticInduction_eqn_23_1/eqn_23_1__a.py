from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_23_1__a(self, B0: float, E: float, r: float, **kwargs):
    # [.pyeqn] E = B0 * a * r ** 2
    result = []
    a = E / (B0 * r**2)
    result.append(a)
    return result
