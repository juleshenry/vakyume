from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_21_4__v(self, B: float, T: float, m: float, q: float, **kwargs):
    # [.pyeqn] q / m = 2 * B * T / (m * v)
    result = []
    v = 2 * B * T / q
    result.append(v)
    return result
