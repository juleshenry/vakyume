from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_7_4ac__n_nc(P_c: float, n_i: float, p: float, p_i: float, **kwargs):
    # [.pyeqn] n_i / n_nc = p_i / (p - P_c)
    result = []
    n_nc = n_i*(-P_c + p)/p_i
    result.append(n_nc)
    return result
