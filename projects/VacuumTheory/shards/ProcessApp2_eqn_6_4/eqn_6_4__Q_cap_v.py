from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_6_4__Q_v(self, delta_h_v: float, w_v: float, **kwargs):
    # [.pyeqn] w_v = 12000 * Q_v / delta_h_v
    result = []
    Q_v = delta_h_v*w_v/12000
    result.append(Q_v)
    return result
