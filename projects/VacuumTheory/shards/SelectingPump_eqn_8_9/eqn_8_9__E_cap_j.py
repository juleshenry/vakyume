from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_8_9__E_j(self, E_m: float, e: float, r: float, s: float, **kwargs):
    # [.pyeqn] r = 2.93 * (E_j * e) / (E_m * s)
    result = []
    E_j = 0.341296928327645*E_m*r*s/e
    result.append(E_j)
    return result
