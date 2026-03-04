from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_2_22__P_s(self, Q_throughput: float, S_p: float, **kwargs):
    # [.pyeqn] Q_throughput = S_p * P_s
    result = []
    P_s = Q_throughput / S_p
    result.append(P_s)
    return result
