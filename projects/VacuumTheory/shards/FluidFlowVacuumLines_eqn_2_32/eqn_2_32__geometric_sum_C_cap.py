from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_2_32__geometric_sum_C(self, C_series: float, **kwargs):
    # [.pyeqn] 1 / C_series = geometric_sum_C
    result = []
    geometric_sum_C = 1/C_series
    result.append(geometric_sum_C)
    return result
