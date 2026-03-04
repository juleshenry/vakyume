from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_4_7__sum_individual_leak_rates(self, W: float, W_T: float, **kwargs):
    # [.pyeqn] W_T = W + sum_individual_leak_rates
    result = []
    sum_individual_leak_rates = -W + W_T
    result.append(sum_individual_leak_rates)
    return result
