from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_3_8__V_P(self, A_C: float, H_2: float, **kwargs):
    # [.pyeqn] V_P = A_C * H_2
    result = []
    V_P = A_C*H_2
    result.append(V_P)
    return result
