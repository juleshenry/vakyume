from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_8_6__K_A(self, K_B: float, W_net: float, **kwargs):
    # [.pyeqn] W_net = K_B - K_A
    result = []
    K_A = K_B - W_net
    result.append(K_A)
    return result
