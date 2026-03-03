from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_3_5__V(self, P: float, P_P: float, V_P: float, **kwargs):
    # [.pyeqn] P_P = P * (V / V_P)
    result = []
    V = P_P*V_P/P
    result.append(V)
    return result
