from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_21_4__q(self, B: float, T: float, m: float, v: float, **kwargs):
    # [.pyeqn] q / m = 2 * B * T / (m * v)
    result = []
    q = 2 * B * T / v
    result.append(q)
    return result
