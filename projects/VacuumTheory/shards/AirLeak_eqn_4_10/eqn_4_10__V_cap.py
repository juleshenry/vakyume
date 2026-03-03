from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_4_10__V(self, T: float, del_P: float, leakage: float, t: float, **kwargs):
    # [.pyeqn] leakage = 0.0059 * V * del_P / t * 530 / T  # lb/hr
    result = []
    V = 0.319795330988168*T*leakage*t/del_P
    result.append(V)
    return result
