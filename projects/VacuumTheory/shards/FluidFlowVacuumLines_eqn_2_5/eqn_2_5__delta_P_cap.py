from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_5__delta_P(D: float, L: float, mu: float, q: float, **kwargs):
    # [.pyeqn] q = 3.141592653589793 * (D ** 4) * delta_P / (128 * L * mu)
    result = []
    delta_P = 40.7436654315252*L*mu*q/D**4
    result.append(delta_P)
    return result
