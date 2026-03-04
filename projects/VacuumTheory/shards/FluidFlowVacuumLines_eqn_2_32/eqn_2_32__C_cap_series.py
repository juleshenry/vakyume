from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_2_32__C_series(self, geometric_sum_C: float, **kwargs):
    # [.pyeqn] 1 / C_series = geometric_sum_C
    result = []
    C_series = 1/geometric_sum_C
    result.append(C_series)
    return result
