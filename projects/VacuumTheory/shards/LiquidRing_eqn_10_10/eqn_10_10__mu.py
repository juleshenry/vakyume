from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_10_10__mu(self, bhp: float, bhp_0: float, rho: float, **kwargs):
    # [.pyeqn] bhp = bhp_0 * (0.5 + 0.0155 * rho ** 0.84 * mu ** 0.16)
    result = []
    mu = (
        -4.7751763343393e-10
        * I
        * (2000.0 * bhp / (bhp_0 * rho**0.84) - 1000.0 / rho**0.84) ** (25 / 4)
    )
    result.append(mu)
    mu = (
        4.7751763343393e-10
        * I
        * (2000.0 * bhp / (bhp_0 * rho**0.84) - 1000.0 / rho**0.84) ** (25 / 4)
    )
    result.append(mu)
    mu = -4.7751763343393e-10 * (
        2000.0 * bhp / (bhp_0 * rho**0.84) - 1000.0 / rho**0.84
    ) ** (25 / 4)
    result.append(mu)
    mu = 4.7751763343393e-10 * (
        2000.0 * bhp / (bhp_0 * rho**0.84) - 1000.0 / rho**0.84
    ) ** (25 / 4)
    result.append(mu)
    return result
