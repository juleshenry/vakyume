from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_5_11__B(self, L_N: float, V_0: float, **kwargs):
    # [.pyeqn] L_N / V_0 = (V_0 + B) / V_0
    result = []
    B = L_N - V_0
    result.append(B)
    return result
