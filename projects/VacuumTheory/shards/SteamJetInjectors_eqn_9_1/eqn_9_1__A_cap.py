from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_9_1__A(rho_s: float, v: float, w_s: float, **kwargs):
    # [.pyeqn] w_s = v * A * rho_s
    result = []
    A = w_s/(rho_s*v)
    result.append(A)
    return result
