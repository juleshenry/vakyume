from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_8_9__v_x(self, L: float, g: float, m: float, x: float, **kwargs):
    # [.pyeqn] L = 0.5 * m * v_x ** 2 + m * g * x
    result = []
    v_x = -sqrt(2.0 * L / m - 2.0 * g * x)
    result.append(v_x)
    v_x = sqrt(2.0 * L / m - 2.0 * g * x)
    result.append(v_x)
    return result
