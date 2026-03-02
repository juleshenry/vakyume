from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_5_10a__V_1(self, D: float, L_0: float, **kwargs):
    # [.pyeqn] L_0 / V_1 = (L_0 / D) / (L_0 / D + 1)
    result = []
    V_1 = D + L_0
    result.append(V_1)
    return result
