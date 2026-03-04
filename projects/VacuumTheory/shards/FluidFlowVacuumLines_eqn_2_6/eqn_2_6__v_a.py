from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_2_6__v_a(self, lambd: float, mu: float, rho: float, **kwargs):
    # [.pyeqn] mu = 0.35 * rho * lambd * v_a
    result = []
    v_a = 2.85714285714286 * mu / (lambd * rho)
    result.append(v_a)
    return result
