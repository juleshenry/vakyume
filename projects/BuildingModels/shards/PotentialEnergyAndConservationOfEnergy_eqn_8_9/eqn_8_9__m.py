from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_8_9__m(self, L: float, g: float, v_x: float, x: float, **kwargs):
    # [.pyeqn] L = 0.5 * m * v_x ** 2 + m * g * x
    result = []
    m = 2.0 * L / (2.0 * g * x + v_x**2)
    result.append(m)
    return result
