from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_5_10a__D(self, L_0: float, V_1: float, **kwargs):
    # [.pyeqn] L_0 / V_1 = (L_0 / D) / (L_0 / D + 1)
    result = []
    D = -L_0 + V_1
    result.append(D)
    return result
