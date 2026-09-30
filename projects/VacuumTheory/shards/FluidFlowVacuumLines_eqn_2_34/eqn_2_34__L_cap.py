from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_2_34__L(
    self, C: float, C_1: float, C_2: float, D: float, P_p: float, mu: float, **kwargs
):
    # [.pyeqn] C = C_1 * (D ** 4 / (mu * L)) * P_p + C_2 * (D ** 3 / L)
    result = []
    L = D**3 * (C_1 * D * P_p + C_2 * mu) / (C * mu)
    result.append(L)
    return result
