from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_5_10c__R(D: float, L_0: float, **kwargs):
    # [.pyeqn] (L_0 / D) / (L_0 / D + 1) = R / (R + 1)
    result = []
    R = L_0/D
    result.append(R)
    return result
