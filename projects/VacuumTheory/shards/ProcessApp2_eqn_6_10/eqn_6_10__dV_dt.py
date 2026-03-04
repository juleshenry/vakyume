from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_6_10__dV_dt(
    self,
    A: float,
    delta_P: float,
    mu: float,
    r_c: float,
    s: float,
    tau: float,
    **kwargs,
):
    # [.pyeqn] dV_dt = (A * delta_P**(1 - s) ) / ( mu * tau * r_c )
    result = []
    dV_dt = A * delta_P ** (1 - s) / (mu * r_c * tau)
    result.append(dV_dt)
    return result
