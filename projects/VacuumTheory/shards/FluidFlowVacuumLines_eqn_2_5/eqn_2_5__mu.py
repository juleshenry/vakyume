from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_2_5__mu(self, D: float, L: float, delta_P: float, q: float, **kwargs):
    # [.pyeqn] q = 3.141592653589793 * (D ** 4) * delta_P / (128 * L * mu)
    result = []
    mu = 0.0245436926061703*D**4*delta_P/(L*q)
    result.append(mu)
    return result
