from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_3_3__P(self, H_1: float, H_2: float, P_P: float, **kwargs):
    # [.pyeqn] P_P - P = H_2 - H_1
    result = []
    P = H_1 - H_2 + P_P
    result.append(P)
    return result
