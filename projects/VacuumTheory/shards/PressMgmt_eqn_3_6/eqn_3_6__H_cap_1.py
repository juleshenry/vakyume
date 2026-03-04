from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_3_6__H_1(self, H_2: float, P: float, V: float, V_P: float, **kwargs):
    # [.pyeqn] P = V_P * (H_2 - H_1) / (V - V_P)
    result = []
    H_1 = H_2 - P*V/V_P + P
    result.append(H_1)
    return result
