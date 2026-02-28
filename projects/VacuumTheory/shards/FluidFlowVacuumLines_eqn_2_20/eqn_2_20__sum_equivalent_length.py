from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_20__sum_equivalent_length(L: float, sum_pipe: float, **kwargs):
    # [.pyeqn] L = sum_pipe + sum_equivalent_length
    result = []
    sum_equivalent_length = L - sum_pipe
    result.append(sum_equivalent_length)
    return result
