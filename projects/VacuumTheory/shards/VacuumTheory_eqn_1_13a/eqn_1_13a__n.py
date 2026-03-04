from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_1_13a__n(self, n_a: float, y_a: float, **kwargs):
    # [.pyeqn] y_a = n_a / n
    result = []
    n = n_a/y_a
    result.append(n)
    return result
