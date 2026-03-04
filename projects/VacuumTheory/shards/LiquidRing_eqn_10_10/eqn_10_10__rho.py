from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_10_10__rho(self, bhp, bhp_0, mu, **kwargs):
    # [.pyeqn] bhp = bhp_0 * (0.5 + 0.0155 * rho ** 0.84 * mu ** 0.16)
    """
    Solves the given equation for rho.
    :param bhp: The brake horsepower (bhp) value.
    :param bhp_0: The brake horsepower at reference conditions.
    :param mu: The friction coefficient.
    :return: The calculated value of rho.
    """
    # Rearrange the equation to isolate rho on one side
    # bhp / bhp_0 = 0.5 + 0.0155 * rho ** 0.84 * mu ** 0.16
    # rho ** 0.84 * mu ** 0.16 = (bhp / bhp_0 - 0.5)
    # rho = ((bhp / bhp_0 - 0.5) / (0.0155 * mu ** 0.16) ** (1 / 0.84)

    # Calculate rho
    rho = ((bhp / bhp_0 - 0.5) ** (1 / 0.84)) / (0.0155 * mu**0.16) ** (1 / 0.84)

    return [rho]
