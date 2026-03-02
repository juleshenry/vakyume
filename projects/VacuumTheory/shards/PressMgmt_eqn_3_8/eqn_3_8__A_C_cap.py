from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_3_8__A_C(self, H_2: float, V_P: float, **kwargs):
    # [.pyeqn] V_P = A_C * H_2
    result = []
    A_C = V_P/H_2
    result.append(A_C)
    return result
