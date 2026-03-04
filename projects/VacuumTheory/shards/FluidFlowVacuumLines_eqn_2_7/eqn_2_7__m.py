from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_2_7__m(self, T: float, k: float, v_a: float, **kwargs):
    # [.pyeqn] v_a = ((8 * k * T) / (3.141592653589793 * m)) ** 0.5
    result = []
    m = 2.54647908947033 * T * k / v_a**2
    result.append(m)
    return result
