from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_4_7__W_T(self, W: float, sum_individual_leak_rates: float, **kwargs):
    # [.pyeqn] W_T = W + sum_individual_leak_rates
    result = []
    W_T = W + sum_individual_leak_rates
    result.append(W_T)
    return result
