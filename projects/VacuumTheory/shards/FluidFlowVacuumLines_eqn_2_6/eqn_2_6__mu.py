from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_2_6__mu(self, lambd: float, rho: float, v_a: float, **kwargs):
    # [.pyeqn] mu = 0.35 * rho * lambd * v_a
    result = []
    mu = 0.35 * lambd * rho * v_a
    result.append(mu)
    return result
