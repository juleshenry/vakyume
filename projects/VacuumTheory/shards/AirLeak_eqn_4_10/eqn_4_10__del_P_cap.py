from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_4_10__del_P(T: float, V: float, leakage: float, t: float, **kwargs):
    # [.pyeqn] leakage = 0.0059 * V * del_P / t * 530 / T  # lb/hr
    result = []
    del_P = 0.319795330988168*T*leakage*t/V
    result.append(del_P)
    return result
