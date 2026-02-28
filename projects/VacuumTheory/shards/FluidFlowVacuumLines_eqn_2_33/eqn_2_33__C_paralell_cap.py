from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_33__C_paralell(arithmetic_sum_C: float, **kwargs):
    # [.pyeqn] 1 / C_paralell = arithmetic_sum_C
    result = []
    C_paralell = 1/arithmetic_sum_C
    result.append(C_paralell)
    return result
