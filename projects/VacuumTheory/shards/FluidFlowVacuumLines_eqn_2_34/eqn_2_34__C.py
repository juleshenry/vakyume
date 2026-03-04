from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_2_34__C(
    self, C_1: float, C_2: float, D: float, L: float, P_p: float, mu: float, **kwargs
):
    # [.pyeqn] C = C_1 * (D ** 4 / (mu * L)) * P_p + C_2 * (D ** 3 / L)
    result = []
    C = D**3 * (C_1 * D * P_p + C_2 * mu) / (L * mu)
    result.append(C)
    return result
