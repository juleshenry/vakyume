from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_4_10__T(V: float, del_P: float, leakage: float, t: float, **kwargs):
    # [.pyeqn] leakage = 0.0059 * V * del_P / t * 530 / T  # lb/hr
    result = []
    T = 3.127*V*del_P/(leakage*t)
    result.append(T)
    return result
