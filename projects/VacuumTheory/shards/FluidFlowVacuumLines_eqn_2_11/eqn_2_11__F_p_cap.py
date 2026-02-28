from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_11__F_p(C_L: float, C_T: float, **kwargs):
    # [.pyeqn] C_T = C_L * F_p
    result = []
    F_p = C_T/C_L
    result.append(F_p)
    return result
