from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_2_34__D(
    self, C: float, C_1: float, C_2: float, L: float, P_p: float, mu: float, **kwargs
):
    # [.pyeqn] C = C_1 * (D ** 4 / (mu * L)) * P_p + C_2 * (D ** 3 / L)
    def _residual(D):
        return (C_1 * (D**4 / (mu * L)) * P_p + C_2 * (D**3 / L)) - (C)

    return [safe_brentq(_residual)]
