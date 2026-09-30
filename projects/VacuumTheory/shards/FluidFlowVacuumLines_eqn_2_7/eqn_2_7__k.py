from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_2_7__k(self, T: float, m: float, v_a: float, **kwargs):
    # [.pyeqn] v_a = ((8 * k * T) / (3.141592653589793 * m)) ** 0.5
    result = []
    k = 0.392699081698724 * m * v_a**2 / T
    result.append(k)
    return result
