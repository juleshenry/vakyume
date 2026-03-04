from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_7_4aa__n_nc(self, n_i: float, p_i: float, p_nc: float, **kwargs):
    # [.pyeqn] n_i / n_nc = p_i / p_nc
    result = []
    n_nc = n_i*p_nc/p_i
    result.append(n_nc)
    return result
