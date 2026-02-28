from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_10_11__T_c(T_s: float, **kwargs):
    # [.pyeqn] T_c = T_s + 10
    result = []
    T_c = T_s + 10
    result.append(T_c)
    return result
