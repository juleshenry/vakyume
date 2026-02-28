from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_1_3__T(k: float, m: float, v: float, **kwargs):
    # [.pyeqn] .5 * m * v**2 = 1.5 * k * T
    result = []
    T = 0.333333333333333*m*v**2/k
    result.append(T)
    return result
