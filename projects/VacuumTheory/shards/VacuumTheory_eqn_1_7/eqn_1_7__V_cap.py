from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_1_7__V(self, R: float, T: float, n: float, p: float, **kwargs):
    # [.pyeqn] p * V = n * R * T
    result = []
    V = R*T*n/p
    result.append(V)
    return result
