from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_8_9__W_C(self, W_NC: float, W_net: float, **kwargs):
    # [.pyeqn] W_net = W_NC + W_C
    result = []
    W_C = -W_NC + W_net
    result.append(W_C)
    return result
