from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_7_15__U(sum_R: float, **kwargs):
    # [.pyeqn] 1 / U = sum_R
    result = []
    U = 1/sum_R
    result.append(U)
    return result
