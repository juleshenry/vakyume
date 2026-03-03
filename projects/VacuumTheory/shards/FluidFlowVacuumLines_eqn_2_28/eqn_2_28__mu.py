from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_2_28__mu(self, C: float, D: float, L: float, P_p: float, **kwargs):
    # [.pyeqn] C = 3.141592653589793 * D ** 4 / (128 * mu * L) * P_p
    result = []
    mu = 0.0245436926061703 * D**4 * P_p / (C * L)
    result.append(mu)
    return result
