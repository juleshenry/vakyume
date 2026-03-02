from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_9_2__rho_s(self, P_m: float, d_n: float, w_s: float, **kwargs):
    # [.pyeqn] w_s = 865.8 * d_n**2 * (P_m * rho_s) ** 0.5
    result = []
    rho_s = 1.334027668054e-6*w_s**2/(P_m*d_n**4)
    result.append(rho_s)
    return result
