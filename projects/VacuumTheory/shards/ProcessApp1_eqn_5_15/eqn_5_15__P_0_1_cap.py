from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_5_15__P_0_1(self, M_1: float, M_2: float, P_0_2: float, a_M_12: float, **kwargs):
    # [.pyeqn] a_M_12 = (P_0_1) / (P_0_2) * (M_2 / M_1) ** 0.4
    result = []
    P_0_1 = P_0_2*a_M_12/(M_2/M_1)**(2/5)
    result.append(P_0_1)
    return result
