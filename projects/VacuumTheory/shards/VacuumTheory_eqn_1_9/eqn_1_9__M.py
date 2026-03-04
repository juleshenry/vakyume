from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_1_9__M(self, P: float, R: float, T: float, rho: float, **kwargs):
    # [.pyeqn] rho = P * M / (R * T)
    result = []
    M = R * T * rho / P
    result.append(M)
    return result
