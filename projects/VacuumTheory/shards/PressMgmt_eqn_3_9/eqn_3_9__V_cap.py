from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_3_9__V(A_C: float, H_1: float, H_2: float, P: float, **kwargs):
    # [.pyeqn] P = A_C * H_2 * (H_2 - H_1) / (V - A_C * H_2)
    result = []
    V = A_C*H_2*(-H_1 + H_2 + P)/P
    result.append(V)
    return result
