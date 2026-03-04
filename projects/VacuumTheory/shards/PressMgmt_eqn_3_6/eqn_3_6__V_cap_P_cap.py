from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_3_6__V_P(self, H_1: float, H_2: float, P: float, V: float, **kwargs):
    # [.pyeqn] P = V_P * (H_2 - H_1) / (V - V_P)
    result = []
    V_P = P*V/(-H_1 + H_2 + P)
    result.append(V_P)
    return result
