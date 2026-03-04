from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_6_9__delta_P(
    self, A: float, dV_dt: float, m: float, mu: float, r: float, r_M: float, **kwargs
):
    # [.pyeqn] dV_dt = (A * delta_P) / (mu * (m / A) * r * delta_P + r_M)
    result = []
    delta_P = A * dV_dt * r_M / (A**2 - dV_dt * m * mu * r)
    result.append(delta_P)
    return result
