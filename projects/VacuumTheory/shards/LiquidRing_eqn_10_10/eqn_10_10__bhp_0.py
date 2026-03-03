from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_10_10__bhp_0(self, bhp: float, mu: float, rho: float, **kwargs):
    # [.pyeqn] bhp = bhp_0 * (0.5 + 0.0155 * rho ** 0.84 * mu ** 0.16)
    result = []
    bhp_0 = 2000.0 * bhp / (31.0 * mu**0.16 * rho**0.84 + 1000.0)
    result.append(bhp_0)
    return result
