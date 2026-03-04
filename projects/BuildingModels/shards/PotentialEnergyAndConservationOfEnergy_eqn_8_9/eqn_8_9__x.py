from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_8_9__x(self, L: float, g: float, m: float, v_x: float, **kwargs):
    # [.pyeqn] L = 0.5 * m * v_x ** 2 + m * g * x
    result = []
    x = (L - 0.5 * m * v_x**2) / (g * m)
    result.append(x)
    return result
