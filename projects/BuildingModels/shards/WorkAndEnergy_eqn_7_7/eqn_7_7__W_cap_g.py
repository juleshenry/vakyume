from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_7_7__W_g(self, W_F: float, W_f: float, W_net: float, **kwargs):
    # [.pyeqn] W_net = W_F + W_g + W_f
    result = []
    W_g = -W_F - W_f + W_net
    result.append(W_g)
    return result
