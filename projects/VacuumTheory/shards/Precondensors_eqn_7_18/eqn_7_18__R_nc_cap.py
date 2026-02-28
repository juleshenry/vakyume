from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_7_18__R_nc(D_0: float, D_LM: float, D_i: float, R_fi: float, R_fo: float, U_0: float, h_c: float, h_i: float, k_w: float, x_w: float, **kwargs):
    # [.pyeqn] 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
    result = []
    R_nc = -D_0*R_fi/D_i - D_0/(D_i*h_i) - D_0*x_w/(D_LM*k_w) - R_fo - 1/h_c + 1/U_0
    result.append(R_nc)
    return result
