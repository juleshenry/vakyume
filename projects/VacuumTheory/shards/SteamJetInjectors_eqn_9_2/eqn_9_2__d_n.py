from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_9_2__d_n(self, P_m: float, rho_s: float, w_s: float, **kwargs):
    # [.pyeqn] w_s = 865.8 * d_n**2 * (P_m * rho_s) ** 0.5
    result = []
    d_n = -0.0339853079285911 * sqrt(w_s / (P_m * rho_s) ** 0.5)
    result.append(d_n)
    d_n = 0.0339853079285911 * sqrt(w_s / (P_m * rho_s) ** 0.5)
    result.append(d_n)
    return result
