from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_2_12__f(
    self, L: float, d: float, delta_P: float, g: float, rho: float, v: float, **kwargs
):
    # [.pyeqn] delta_P = 4.31 * rho * f * L * v ** 2 / (2 * d * g)
    result = []
    f = 0.464037122969838 * d * delta_P * g / (L * rho * v**2)
    result.append(f)
    return result
