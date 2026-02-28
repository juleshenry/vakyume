from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_1_7__T(R: float, V: float, n: float, p: float, **kwargs):
    # [.pyeqn] p * V = n * R * T
    result = []
    T = V*p/(R*n)
    result.append(T)
    return result
