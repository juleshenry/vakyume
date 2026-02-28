from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_32__geometric_sum_C(C_series: float, **kwargs):
    # [.pyeqn] 1 / C_series = geometric_sum_C
    result = []
    geometric_sum_C = 1/C_series
    result.append(geometric_sum_C)
    return result
