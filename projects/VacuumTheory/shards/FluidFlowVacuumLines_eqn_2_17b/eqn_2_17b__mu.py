from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_2_17b__mu(self, L: float, d: float, delta_P: float, q: float, **kwargs):
    # [.pyeqn] delta_P = 0.105 * mu * L * q / d**4
    result = []
    mu = 9.52380952380952 * d**4 * delta_P / (L * q)
    result.append(mu)
    return result
