from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_3_11__P(self, A_C: float, H_2: float, V: float, **kwargs):
    # [.pyeqn] P = A_C / V * (H_2) ** 2
    result = []
    P = A_C*H_2**2/V
    result.append(P)
    return result
