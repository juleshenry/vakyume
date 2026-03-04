from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_7_1__y_i(self, P: float, p_i: float, **kwargs):
    # [.pyeqn] y_i = p_i / P
    result = []
    y_i = p_i/P
    result.append(y_i)
    return result
