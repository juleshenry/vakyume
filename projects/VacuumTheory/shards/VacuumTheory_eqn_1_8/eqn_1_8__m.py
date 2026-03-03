from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_1_8__m(self, M: float, P: float, R: float, T: float, V: float, **kwargs):
    # [.pyeqn] P * V = m / M * R * T
    result = []
    m = M*P*V/(R*T)
    result.append(m)
    return result
