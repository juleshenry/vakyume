from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_3_6__P(self, H_1: float, H_2: float, V: float, V_P: float, **kwargs):
    # [.pyeqn] P = V_P * (H_2 - H_1) / (V - V_P)
    result = []
    P = V_P*(-H_1 + H_2)/(V - V_P)
    result.append(P)
    return result
