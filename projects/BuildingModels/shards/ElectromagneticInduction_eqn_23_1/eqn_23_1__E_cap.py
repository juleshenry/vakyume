from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_23_1__E(self, B0: float, a: float, r: float, **kwargs):
    # [.pyeqn] E = B0 * a * r ** 2
    result = []
    E = B0 * a * r**2
    result.append(E)
    return result
