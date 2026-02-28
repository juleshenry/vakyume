from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_1_8__m(M: float, P: float, R: float, T: float, V: float, **kwargs):
    # [.pyeqn] P * V = m / M * R * T
    result = []
    m = M*P*V/(R*T)
    result.append(m)
    return result
