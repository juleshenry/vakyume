from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_7_7__W_F(self, W_f: float, W_g: float, W_net: float, **kwargs):
    # [.pyeqn] W_net = W_F + W_g + W_f
    result = []
    W_F = -W_f - W_g + W_net
    result.append(W_F)
    return result
