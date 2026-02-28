from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_3_5__V(P: float, P_P: float, V_P: float, **kwargs):
    # [.pyeqn] P_P = P * (V / V_P)
    result = []
    V = P_P*V_P/P
    result.append(V)
    return result
