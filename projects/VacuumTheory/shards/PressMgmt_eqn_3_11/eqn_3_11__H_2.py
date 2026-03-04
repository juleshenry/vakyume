from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_3_11__H_2(self, A_C: float, P: float, V: float, **kwargs):
    # [.pyeqn] P = A_C / V * (H_2) ** 2
    result = []
    H_2 = -sqrt(P*V/A_C)
    result.append(H_2)
    H_2 = sqrt(P*V/A_C)
    result.append(H_2)
    return result
