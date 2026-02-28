from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_11_3__F_s(t: float, t_c: float, **kwargs):
    # [.pyeqn] t = t_c * F_s
    result = []
    F_s = t/t_c
    result.append(F_s)
    return result
