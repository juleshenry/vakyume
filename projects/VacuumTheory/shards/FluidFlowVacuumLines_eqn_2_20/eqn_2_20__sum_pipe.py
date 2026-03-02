from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_20__sum_pipe(self, L: float, sum_equivalent_length: float, **kwargs):
    # [.pyeqn] L = sum_pipe + sum_equivalent_length
    result = []
    sum_pipe = L - sum_equivalent_length
    result.append(sum_pipe)
    return result
