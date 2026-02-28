from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_7_7__M(P: float, P_i_0: float, W_air: float, W_i: float, epsilon_i: float, p_c: float, x_i: float, **kwargs):
    # [.pyeqn] W_i = W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c))
    result = []
    M = 29*W_i*(P - p_c)/(P_i_0*W_air*epsilon_i*x_i)
    result.append(M)
    return result
