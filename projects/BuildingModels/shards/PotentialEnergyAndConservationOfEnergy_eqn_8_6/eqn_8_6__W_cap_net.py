from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_8_6__W_net(self, K_A: float, K_B: float, **kwargs):
    # [.pyeqn] W_net = K_B - K_A
    result = []
    W_net = -K_A + K_B
    result.append(W_net)
    return result
