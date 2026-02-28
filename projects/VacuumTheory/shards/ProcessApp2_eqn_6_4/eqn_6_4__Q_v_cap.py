from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_6_4__Q_v(delta_h_v: float, w_v: float, **kwargs):
    # [.pyeqn] w_v = 12000 * Q_v / delta_h_v
    result = []
    Q_v = delta_h_v*w_v/12000
    result.append(Q_v)
    return result
