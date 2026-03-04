from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_5_4__P(self, P_0_i: float, x_i: float, y_i: float, **kwargs):
    # [.pyeqn] y_i * P = x_i * P_0_i
    result = []
    P = P_0_i*x_i/y_i
    result.append(P)
    return result
