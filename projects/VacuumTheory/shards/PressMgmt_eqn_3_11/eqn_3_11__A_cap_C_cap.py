from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_3_11__A_C(self, H_2: float, P: float, V: float, **kwargs):
    # [.pyeqn] P = A_C / V * (H_2) ** 2
    result = []
    A_C = P*V/H_2**2
    result.append(A_C)
    return result
