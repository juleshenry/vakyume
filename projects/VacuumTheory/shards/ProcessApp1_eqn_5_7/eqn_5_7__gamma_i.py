from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_5_7__gamma_i(self, P: float, P_0_i: float, x_i: float, y_i: float, **kwargs):
    # [.pyeqn] y_i * P = x_i * gamma_i * P_0_i
    result = []
    gamma_i = P*y_i/(P_0_i*x_i)
    result.append(gamma_i)
    return result
