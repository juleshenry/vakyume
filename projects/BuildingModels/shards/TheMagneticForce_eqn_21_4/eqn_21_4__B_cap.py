from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_21_4__B(self, T: float, m: float, q: float, v: float, **kwargs):
    # [.pyeqn] q / m = 2 * B * T / (m * v)
    result = []
    B = q * v / (2 * T)
    result.append(B)
    return result
