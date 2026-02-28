from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_3_6__V_P(H_1: float, H_2: float, P: float, V: float, **kwargs):
    # [.pyeqn] P = V_P * (H_2 - H_1) / (V - V_P)
    result = []
    V_P = P*V/(-H_1 + H_2 + P)
    result.append(V_P)
    return result
