from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_3_6__V(self, H_1: float, H_2: float, P: float, V_P: float, **kwargs):
    # [.pyeqn] P = V_P * (H_2 - H_1) / (V - V_P)
    result = []
    V = V_P*(-H_1 + H_2 + P)/P
    result.append(V)
    return result
