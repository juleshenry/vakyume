from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_26__P_p(self, D: float, L: float, P_downstream: float, P_upstream: float, mu: float, q: float, **kwargs):
    # [.pyeqn] q * P_p = 3.141592653589793 * D ** 4 / (128 * mu * L) * P_p * (P_upstream - P_downstream)
    result = []
    P_p = 0.0
    result.append(P_p)
    return result
