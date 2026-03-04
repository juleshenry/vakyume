from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_7_4aa__p_nc(self, n_i: float, n_nc: float, p_i: float, **kwargs):
    # [.pyeqn] n_i / n_nc = p_i / p_nc
    result = []
    p_nc = n_nc*p_i/n_i
    result.append(p_nc)
    return result
