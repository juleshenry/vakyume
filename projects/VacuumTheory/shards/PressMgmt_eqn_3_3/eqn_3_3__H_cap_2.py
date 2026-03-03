from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_3_3__H_2(self, H_1: float, P: float, P_P: float, **kwargs):
    # [.pyeqn] P_P - P = H_2 - H_1
    result = []
    H_2 = H_1 - P + P_P
    result.append(H_2)
    return result
