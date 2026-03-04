from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_17_3__E(self, Q: float, R: float, **kwargs):
    # [.pyeqn] E = Q / (4 * R ** 2)
    result = []
    E = Q / (4 * R**2)
    result.append(E)
    return result
