from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_2_28__D(self, C: float, L: float, P_p: float, mu: float, **kwargs):
    # [.pyeqn] C = 3.141592653589793 * D ** 4 / (128 * mu * L) * P_p
    result = []
    D = -2.52647511098426 * I * (C * L * mu / P_p) ** (1 / 4)
    result.append(D)
    D = 2.52647511098426 * I * (C * L * mu / P_p) ** (1 / 4)
    result.append(D)
    D = -2.52647511098426 * (C * L * mu / P_p) ** (1 / 4)
    result.append(D)
    D = 2.52647511098426 * (C * L * mu / P_p) ** (1 / 4)
    result.append(D)
    return result
