from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_10_21__P_prime(self, P: float, P_d: float, **kwargs):
    # [.pyeqn] P_prime = P / P_d * 760
    result = []
    P_prime = 760*P/P_d
    result.append(P_prime)
    return result
