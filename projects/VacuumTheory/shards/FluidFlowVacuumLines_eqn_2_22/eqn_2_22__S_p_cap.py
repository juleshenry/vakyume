from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_22__S_p(P_s: float, Q_throughput: float, **kwargs):
    # [.pyeqn] Q_throughput = S_p * P_s
    result = []
    S_p = Q_throughput/P_s
    result.append(S_p)
    return result
