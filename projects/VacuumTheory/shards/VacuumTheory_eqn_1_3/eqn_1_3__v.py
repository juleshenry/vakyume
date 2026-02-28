from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_1_3__v(T: float, k: float, m: float, **kwargs):
    # [.pyeqn] .5 * m * v**2 = 1.5 * k * T
    result = []
    v = -1.73205080756888*sqrt(T*k/m)
    result.append(v)
    v = 1.73205080756888*sqrt(T*k/m)
    result.append(v)
    return result
