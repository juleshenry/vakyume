from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_6_4__delta_h_v(self, Q_v: float, w_v: float, **kwargs):
    # [.pyeqn] w_v = 12000 * Q_v / delta_h_v
    result = []
    delta_h_v = 12000*Q_v/w_v
    result.append(delta_h_v)
    return result
