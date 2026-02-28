from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_7_4ac__p(P_c: float, n_i: float, n_nc: float, p_i: float, **kwargs):
    # [.pyeqn] n_i / n_nc = p_i / (p - P_c)
    result = []
    p = P_c + n_nc*p_i/n_i
    result.append(p)
    return result
