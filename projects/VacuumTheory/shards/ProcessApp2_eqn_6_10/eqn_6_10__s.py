from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_6_10__s(
    self,
    A: float,
    dV_dt: float,
    delta_P: float,
    mu: float,
    r_c: float,
    tau: float,
    **kwargs,
):
    # [.pyeqn] dV_dt = (A * delta_P**(1 - s) ) / ( mu * tau * r_c )
    result = []
    s = log(A * delta_P / (dV_dt * mu * r_c * tau)) / log(delta_P)
    result.append(s)
    return result
