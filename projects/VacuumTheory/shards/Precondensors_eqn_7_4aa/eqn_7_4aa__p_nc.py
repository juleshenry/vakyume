from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_7_4aa__p_nc(n_i: float, n_nc: float, p_i: float, **kwargs):
    # [.pyeqn] n_i / n_nc = p_i / p_nc
    result = []
    p_nc = n_nc*p_i/n_i
    result.append(p_nc)
    return result
