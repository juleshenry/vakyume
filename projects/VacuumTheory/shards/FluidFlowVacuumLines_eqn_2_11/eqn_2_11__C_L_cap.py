from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_11__C_L(self, C_T: float, F_p: float, **kwargs):
    # [.pyeqn] C_T = C_L * F_p
    result = []
    C_L = C_T/F_p
    result.append(C_L)
    return result
