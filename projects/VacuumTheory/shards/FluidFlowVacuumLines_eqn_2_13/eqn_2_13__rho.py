from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_13__rho(L: float, d: float, delta_P: float, f: float, q: float, **kwargs):
    # [.pyeqn] delta_P = 2.15 * rho * f * L * q ** 2 / (d ** 5)
    result = []
    rho = 0.465116279069767*d**5*delta_P/(L*f*q**2)
    result.append(rho)
    return result
