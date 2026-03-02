from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_25__P_1(self, C: float, P_2: float, Q_throughput: float, **kwargs):
    # [.pyeqn] C = Q_throughput / (P_1 - P_2)
    result = []
    P_1 = P_2 + Q_throughput/C
    result.append(P_1)
    return result
