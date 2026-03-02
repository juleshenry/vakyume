from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_4_7__W_T(self, W: float, sum_individual_leak_rates: float, **kwargs):
    # [.pyeqn] W_T = W + sum_individual_leak_rates
    result = []
    W_T = W + sum_individual_leak_rates
    result.append(W_T)
    return result
