from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_3_5__V_P(self, P: float, P_P: float, V: float, **kwargs):
    # [.pyeqn] P_P = P * (V / V_P)
    result = []
    V_P = P*V/P_P
    result.append(V_P)
    return result
