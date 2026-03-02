from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_12__delta_P(self, L: float, d: float, f: float, g: float, rho: float, v: float, **kwargs):
    # [.pyeqn] delta_P = 4.31 * rho * f * L * v ** 2 / (2 * d * g)
    result = []
    delta_P = 2.155*L*f*rho*v**2/(d*g)
    result.append(delta_P)
    return result
