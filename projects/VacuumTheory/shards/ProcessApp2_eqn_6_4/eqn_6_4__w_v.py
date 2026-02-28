from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_6_4__w_v(Q_v: float, delta_h_v: float, **kwargs):
    # [.pyeqn] w_v = 12000 * Q_v / delta_h_v
    result = []
    w_v = 12000*Q_v/delta_h_v
    result.append(w_v)
    return result
