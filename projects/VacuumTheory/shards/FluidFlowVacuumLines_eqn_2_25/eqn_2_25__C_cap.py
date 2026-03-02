from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_25__C(self, P_1: float, P_2: float, Q_throughput: float, **kwargs):
    # [.pyeqn] C = Q_throughput / (P_1 - P_2)
    result = []
    C = Q_throughput/(P_1 - P_2)
    result.append(C)
    return result
