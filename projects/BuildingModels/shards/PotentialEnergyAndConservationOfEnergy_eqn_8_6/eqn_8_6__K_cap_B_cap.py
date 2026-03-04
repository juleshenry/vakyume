from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_8_6__K_B(self, K_A: float, W_net: float, **kwargs):
    # [.pyeqn] W_net = K_B - K_A
    result = []
    K_B = K_A + W_net
    result.append(K_B)
    return result
