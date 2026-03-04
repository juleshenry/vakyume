from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_1_3__T(self, k: float, m: float, v: float, **kwargs):
    # [.pyeqn] .5 * m * v**2 = 1.5 * k * T
    result = []
    T = 0.333333333333333*m*v**2/k
    result.append(T)
    return result
