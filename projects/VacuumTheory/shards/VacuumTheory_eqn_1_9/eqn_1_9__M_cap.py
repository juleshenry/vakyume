from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_1_9__M(P: float, R: float, T: float, rho: float, **kwargs):
    # [.pyeqn] rho = P * M / (R * T)
    result = []
    M = R*T*rho/P
    result.append(M)
    return result
