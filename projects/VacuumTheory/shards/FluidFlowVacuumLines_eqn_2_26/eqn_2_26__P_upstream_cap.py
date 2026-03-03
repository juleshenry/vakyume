from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_26__P_upstream(self, D: float, L: float, P_downstream: float, P_p: float, mu: float, q: float, **kwargs):
    # [.pyeqn] q * P_p = 3.141592653589793 * D ** 4 / (128 * mu * L) * P_p * (P_upstream - P_downstream)
    result = []
    P_upstream = P_downstream + 40.7436654315252*L*mu*q/D**4
    result.append(P_upstream)
    return result
