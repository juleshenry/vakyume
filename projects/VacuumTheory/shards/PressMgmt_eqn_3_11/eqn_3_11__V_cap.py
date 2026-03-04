from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_3_11__V(self, A_C: float, H_2: float, P: float, **kwargs):
    # [.pyeqn] P = A_C / V * (H_2) ** 2
    result = []
    V = A_C*H_2**2/P
    result.append(V)
    return result
