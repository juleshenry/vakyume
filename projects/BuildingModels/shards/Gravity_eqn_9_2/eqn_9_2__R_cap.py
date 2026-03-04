from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_9_2__R(self, G: float, M: float, T: float, pi: float, **kwargs):
    # [.pyeqn] T ** 2 = G * M * R ** 3 / (4 * pi ** 2)
    result = []
    R = 2 ** (2 / 3) * (T**2 * pi**2 / (G * M)) ** (1 / 3)
    result.append(R)
    R = (
        -(2 ** (2 / 3)) * (T**2 * pi**2 / (G * M)) ** (1 / 3) / 2
        - 2 ** (2 / 3) * sqrt(3) * I * (T**2 * pi**2 / (G * M)) ** (1 / 3) / 2
    )
    result.append(R)
    R = (
        -(2 ** (2 / 3)) * (T**2 * pi**2 / (G * M)) ** (1 / 3) / 2
        + 2 ** (2 / 3) * sqrt(3) * I * (T**2 * pi**2 / (G * M)) ** (1 / 3) / 2
    )
    result.append(R)
    return result
