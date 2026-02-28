from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_8_9__s(E_j: float, E_m: float, e: float, r: float, **kwargs):
    # [.pyeqn] r = 2.93 * (E_j * e) / (E_m * s)
    result = []
    s = 2.93*E_j*e/(E_m*r)
    result.append(s)
    return result
