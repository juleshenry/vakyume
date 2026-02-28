from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_25__Q_throughput(C: float, P_1: float, P_2: float, **kwargs):
    # [.pyeqn] C = Q_throughput / (P_1 - P_2)
    result = []
    Q_throughput = C*(P_1 - P_2)
    result.append(Q_throughput)
    return result
