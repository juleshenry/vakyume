from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_3_9__P(self, A_C: float, H_1: float, H_2: float, V: float, **kwargs):
    # [.pyeqn] P = A_C * H_2 * (H_2 - H_1) / (V - A_C * H_2)
    result = []
    P = A_C*H_2*(H_1 - H_2)/(A_C*H_2 - V)
    result.append(P)
    return result
