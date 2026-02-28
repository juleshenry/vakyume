from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_5__q(D: float, L: float, delta_P: float, mu: float, **kwargs):
    # [.pyeqn] q = 3.141592653589793 * (D ** 4) * delta_P / (128 * L * mu)
    result = []
    q = 0.0245436926061703*D**4*delta_P/(L*mu)
    result.append(q)
    return result
