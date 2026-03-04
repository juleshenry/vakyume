from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_2_26__mu(self, D: float, L: float, P_downstream: float, P_p: float, P_upstream: float, q: float, **kwargs):
    # [.pyeqn] q * P_p = 3.141592653589793 * D ** 4 / (128 * mu * L) * P_p * (P_upstream - P_downstream)
    result = []
    mu = 0.0245436926061703*D**4*(-P_downstream + P_upstream)/(L*q)
    result.append(mu)
    return result
