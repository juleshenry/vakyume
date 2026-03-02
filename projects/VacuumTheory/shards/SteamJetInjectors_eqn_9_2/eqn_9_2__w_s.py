from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_9_2__w_s(self, P_m: float, d_n: float, rho_s: float, **kwargs):
    # [.pyeqn] w_s = 865.8 * d_n**2 * (P_m * rho_s) ** 0.5
    result = []
    w_s = 865.8*d_n**2*sqrt(P_m*rho_s)
    result.append(w_s)
    return result
