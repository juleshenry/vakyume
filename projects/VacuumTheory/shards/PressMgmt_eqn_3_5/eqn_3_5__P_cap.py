from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_3_5__P(P_P: float, V: float, V_P: float, **kwargs):
    # [.pyeqn] P_P = P * (V / V_P)
    result = []
    P = P_P*V_P/V
    result.append(P)
    return result
