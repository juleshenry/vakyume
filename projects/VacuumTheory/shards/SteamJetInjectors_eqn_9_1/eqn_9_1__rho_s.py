from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_9_1__rho_s(self, A: float, v: float, w_s: float, **kwargs):
    # [.pyeqn] w_s = v * A * rho_s
    result = []
    rho_s = w_s/(A*v)
    result.append(rho_s)
    return result
