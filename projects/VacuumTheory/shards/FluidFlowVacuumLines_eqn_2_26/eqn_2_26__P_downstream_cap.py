from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_26__P_downstream(D: float, L: float, P_p: float, P_upstream: float, mu: float, q: float, **kwargs):
    # [.pyeqn] q * P_p = 3.141592653589793 * D ** 4 / (128 * mu * L) * P_p * (P_upstream - P_downstream)
    result = []
    P_downstream = P_upstream - 40.7436654315252*L*mu*q/D**4
    result.append(P_downstream)
    return result
