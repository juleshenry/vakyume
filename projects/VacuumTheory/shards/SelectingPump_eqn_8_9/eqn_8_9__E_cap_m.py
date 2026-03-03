from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_8_9__E_m(self, E_j: float, e: float, r: float, s: float, **kwargs):
    # [.pyeqn] r = 2.93 * (E_j * e) / (E_m * s)
    result = []
    E_m = 2.93*E_j*e/(r*s)
    result.append(E_m)
    return result
