from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_7_1__p_i(self, P: float, y_i: float, **kwargs):
    # [.pyeqn] y_i = p_i / P
    result = []
    p_i = P*y_i
    result.append(p_i)
    return result
