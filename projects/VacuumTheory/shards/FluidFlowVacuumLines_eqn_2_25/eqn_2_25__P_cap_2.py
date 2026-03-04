from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_2_25__P_2(self, C: float, P_1: float, Q_throughput: float, **kwargs):
    # [.pyeqn] C = Q_throughput / (P_1 - P_2)
    result = []
    P_2 = P_1 - Q_throughput / C
    result.append(P_2)
    return result
