from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_33__arithmetic_sum_C(C_paralell: float, **kwargs):
    # [.pyeqn] 1 / C_paralell = arithmetic_sum_C
    result = []
    arithmetic_sum_C = 1/C_paralell
    result.append(arithmetic_sum_C)
    return result
