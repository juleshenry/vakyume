from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_4_7__W(self, W_T: float, sum_individual_leak_rates: float, **kwargs):
    # [.pyeqn] W_T = W + sum_individual_leak_rates
    result = []
    W = W_T - sum_individual_leak_rates
    result.append(W)
    return result
