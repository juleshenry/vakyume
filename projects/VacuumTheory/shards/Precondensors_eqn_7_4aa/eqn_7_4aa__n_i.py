from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_7_4aa__n_i(n_nc: float, p_i: float, p_nc: float, **kwargs):
    # [.pyeqn] n_i / n_nc = p_i / p_nc
    result = []
    n_i = n_nc*p_i/p_nc
    result.append(n_i)
    return result
