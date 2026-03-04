from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_1_7__T(self, R: float, V: float, n: float, p: float, **kwargs):
    # [.pyeqn] p * V = n * R * T
    result = []
    T = V*p/(R*n)
    result.append(T)
    return result
