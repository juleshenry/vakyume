from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_28__C(self, D: float, L: float, P_p: float, mu: float, **kwargs):
    # [.pyeqn] C = 3.141592653589793 * D ** 4 / (128 * mu * L) * P_p
    result = []
    C = 0.0245436926061703*D**4*P_p/(L*mu)
    result.append(C)
    return result
