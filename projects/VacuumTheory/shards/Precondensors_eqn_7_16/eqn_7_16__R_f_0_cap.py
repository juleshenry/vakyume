from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_7_16__R_f_0(D_0: float, D_LM: float, D_i: float, R_fi: float, U_0: float, h_0: float, h_i: float, k_w: float, x_w: float, **kwargs):
    # [.pyeqn] 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
    result = []
    R_f_0 = -D_0*R_fi/D_i - D_0/(D_i*h_i) - D_0*x_w/(D_LM*k_w) - 1/h_0 + 1/U_0
    result.append(R_f_0)
    return result
