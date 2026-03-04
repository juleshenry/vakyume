from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_2_17b__q(self, L: float, d: float, delta_P: float, mu: float, **kwargs):
    # [.pyeqn] delta_P = 0.105 * mu * L * q / d**4
    result = []
    q = 9.52380952380952 * d**4 * delta_P / (L * mu)
    result.append(q)
    return result
