from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_9_1__v(self, A: float, rho_s: float, w_s: float, **kwargs):
    # [.pyeqn] w_s = v * A * rho_s
    result = []
    v = w_s/(A*rho_s)
    result.append(v)
    return result
