from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_5_2a__alpha_1_2(self, K_1: float, K_2: float, **kwargs):
    # [.pyeqn] alpha_1_2 = K_1 / K_2
    result = []
    alpha_1_2 = K_1 / K_2
    result.append(alpha_1_2)
    return result
