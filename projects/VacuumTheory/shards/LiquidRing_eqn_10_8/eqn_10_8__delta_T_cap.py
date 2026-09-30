from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_10_8__delta_T(
    self,
    bhp: float,
    c_p: float,
    delta_h_i: float,
    f_a: float,
    rho: float,
    w_i: float,
    **kwargs,
):
    # [.pyeqn] delta_T = (2545 * bhp + w_i * delta_h_i) / ( 8.02 * f_a * rho * c_p )
    result = []
    delta_T = 0.124688279301746 * (2545.0 * bhp + delta_h_i * w_i) / (c_p * f_a * rho)
    result.append(delta_T)
    return result
