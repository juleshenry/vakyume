from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_1_8__P(self, M: float, R: float, T: float, V: float, m: float, **kwargs):
    # [.pyeqn] P * V = m / M * R * T
    result = []
    P = R*T*m/(M*V)
    result.append(P)
    return result
