from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_1_7__R(self, T: float, V: float, n: float, p: float, **kwargs):
    # [.pyeqn] p * V = n * R * T
    result = []
    R = V*p/(T*n)
    result.append(R)
    return result
