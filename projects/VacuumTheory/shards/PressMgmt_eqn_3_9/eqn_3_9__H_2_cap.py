from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_3_9__H_2(A_C: float, H_1: float, P: float, V: float, **kwargs):
    # [.pyeqn] P = A_C * H_2 * (H_2 - H_1) / (V - A_C * H_2)
    result = []
    H_2 = (A_C*(H_1 - P) - sqrt(A_C*(A_C*H_1**2 - 2*A_C*H_1*P + A_C*P**2 + 4*P*V)))/(2*A_C)
    result.append(H_2)
    H_2 = (A_C*(H_1 - P) + sqrt(A_C*(A_C*H_1**2 - 2*A_C*H_1*P + A_C*P**2 + 4*P*V)))/(2*A_C)
    result.append(H_2)
    return result
