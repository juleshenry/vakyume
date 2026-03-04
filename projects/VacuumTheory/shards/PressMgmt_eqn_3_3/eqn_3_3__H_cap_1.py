from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_3_3__H_1(self, H_2: float, P: float, P_P: float, **kwargs):
    # [.pyeqn] P_P - P = H_2 - H_1
    result = []
    H_1 = H_2 + P - P_P
    result.append(H_1)
    return result
