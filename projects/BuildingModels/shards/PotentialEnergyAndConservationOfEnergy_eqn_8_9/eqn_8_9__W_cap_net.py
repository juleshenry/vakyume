from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_8_9__W_net(self, W_C: float, W_NC: float, **kwargs):
    # [.pyeqn] W_net = W_NC + W_C
    result = []
    W_net = W_C + W_NC
    result.append(W_net)
    return result
