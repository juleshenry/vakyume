from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_8_9__g(self, L: float, m: float, v_x: float, x: float, **kwargs):
    # [.pyeqn] L = 0.5 * m * v_x ** 2 + m * g * x
    result = []
    g = (L - 0.5 * m * v_x**2) / (m * x)
    result.append(g)
    return result
