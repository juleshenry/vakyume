from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_1_9__P(self, M: float, R: float, T: float, rho: float, **kwargs):
    # [.pyeqn] rho = P * M / (R * T)
    result = []
    P = R*T*rho/M
    result.append(P)
    return result
