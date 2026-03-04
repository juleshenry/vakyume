from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_10_21__P_d(self, P: float, P_prime: float, **kwargs):
    # [.pyeqn] P_prime = P / P_d * 760
    result = []
    P_d = 760*P/P_prime
    result.append(P_d)
    return result
