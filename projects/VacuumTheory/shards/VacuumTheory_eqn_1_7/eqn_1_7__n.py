from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_1_7__n(self, R: float, T: float, V: float, p: float, **kwargs):
    # [.pyeqn] p * V = n * R * T
    result = []
    n = V*p/(R*T)
    result.append(n)
    return result
