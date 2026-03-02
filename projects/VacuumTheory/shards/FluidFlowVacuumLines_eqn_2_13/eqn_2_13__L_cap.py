from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_13__L(self, d: float, delta_P: float, f: float, q: float, rho: float, **kwargs):
    # [.pyeqn] delta_P = 2.15 * rho * f * L * q ** 2 / (d ** 5)
    result = []
    L = 0.465116279069767*d**5*delta_P/(f*q**2*rho)
    result.append(L)
    return result
