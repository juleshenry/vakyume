from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_3_9__H_1(self, A_C: float, H_2: float, P: float, V: float, **kwargs):
    # [.pyeqn] P = A_C * H_2 * (H_2 - H_1) / (V - A_C * H_2)
    result = []
    H_1 = H_2 + P - P*V/(A_C*H_2)
    result.append(H_1)
    return result
