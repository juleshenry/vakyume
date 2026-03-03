from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_2_34__P_p(
    self, C: float, C_1: float, C_2: float, D: float, L: float, mu: float, **kwargs
):
    # [.pyeqn] C = C_1 * (D ** 4 / (mu * L)) * P_p + C_2 * (D ** 3 / L)
    result = []
    P_p = mu * (C * L - C_2 * D**3) / (C_1 * D**4)
    result.append(P_p)
    return result
