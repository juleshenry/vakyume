from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_3_3__H_1(self, H_2: float, P: float, P_P: float, **kwargs):
    # [.pyeqn] P_P - P = H_2 - H_1
    result = []
    H_1 = H_2 + P - P_P
    result.append(H_1)
    return result
