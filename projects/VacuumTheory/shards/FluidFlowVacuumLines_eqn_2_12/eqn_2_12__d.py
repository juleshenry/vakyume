from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_12__d(self, L: float, delta_P: float, f: float, g: float, rho: float, v: float, **kwargs):
    # [.pyeqn] delta_P = 4.31 * rho * f * L * v ** 2 / (2 * d * g)
    result = []
    d = 2.155*L*f*rho*v**2/(delta_P*g)
    result.append(d)
    return result
