from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_3_8__H_2(A_C: float, V_P: float, **kwargs):
    # [.pyeqn] V_P = A_C * H_2
    result = []
    H_2 = V_P/A_C
    result.append(H_2)
    return result
