from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_7_1__P(self, p_i: float, y_i: float, **kwargs):
    # [.pyeqn] y_i = p_i / P
    result = []
    P = p_i/y_i
    result.append(P)
    return result
