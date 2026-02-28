from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_1_8__V(M: float, P: float, R: float, T: float, m: float, **kwargs):
    # [.pyeqn] P * V = m / M * R * T
    result = []
    V = R*T*m/(M*P)
    result.append(V)
    return result
