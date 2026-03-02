from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_5_11__V_0(self, B: float, L_N: float, **kwargs):
    # [.pyeqn] L_N / V_0 = (V_0 + B) / V_0
    result = []
    V_0 = -B + L_N
    result.append(V_0)
    return result
