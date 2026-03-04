from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_2_22__S_p(self, P_s: float, Q_throughput: float, **kwargs):
    # [.pyeqn] Q_throughput = S_p * P_s
    result = []
    S_p = Q_throughput / P_s
    result.append(S_p)
    return result
