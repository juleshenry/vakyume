from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_2_22__Q_throughput(self, P_s: float, S_p: float, **kwargs):
    # [.pyeqn] Q_throughput = S_p * P_s
    result = []
    Q_throughput = P_s * S_p
    result.append(Q_throughput)
    return result
