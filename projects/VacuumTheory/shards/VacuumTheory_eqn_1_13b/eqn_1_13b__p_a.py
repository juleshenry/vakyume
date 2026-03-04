from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_1_13b__p_a(self, P: float, y_a: float, **kwargs):
    # [.pyeqn] y_a = p_a / P
    result = []
    p_a = P * y_a
    result.append(p_a)
    return result
