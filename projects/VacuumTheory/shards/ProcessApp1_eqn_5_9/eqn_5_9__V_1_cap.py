from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_5_9__V_1(self, D: float, L_0: float, **kwargs):
    # [.pyeqn] L_0 / V_1 = L_0 / (L_0 + D)
    result = []
    V_1 = D + L_0
    result.append(V_1)
    return result
