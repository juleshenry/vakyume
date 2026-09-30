from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_5_14__P_0(self, M: float, T: float, W_E: float, **kwargs):
    # [.pyeqn] W_E = 0.0583 * P_0 * (M / T) ** 0.5
    result = []
    P_0 = 17.1526586620926 * W_E / sqrt(M / T)
    result.append(P_0)
    return result
