from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_2_17__d(self, L: float, delta_P: float, mu: float, v: float, **kwargs):
    # [.pyeqn] delta_P = 0.0345* mu * L * v / d**2
    result = []
    d = -0.185741756210067 * sqrt(L * mu * v / delta_P)
    result.append(d)
    d = 0.185741756210067 * sqrt(L * mu * v / delta_P)
    result.append(d)
    return result
