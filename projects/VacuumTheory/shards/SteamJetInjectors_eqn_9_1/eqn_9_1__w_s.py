from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_9_1__w_s(A: float, rho_s: float, v: float, **kwargs):
    # [.pyeqn] w_s = v * A * rho_s
    result = []
    w_s = A*rho_s*v
    result.append(w_s)
    return result
