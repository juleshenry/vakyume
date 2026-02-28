from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_28__P_p(C: float, D: float, L: float, mu: float, **kwargs):
    # [.pyeqn] C = 3.141592653589793 * D ** 4 / (128 * mu * L) * P_p
    result = []
    P_p = 40.7436654315252*C*L*mu/D**4
    result.append(P_p)
    return result
