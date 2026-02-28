from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_10_12__T_c(T_s: float, **kwargs):
    # [.pyeqn] T_c = T_s + 5
    result = []
    T_c = T_s + 5
    result.append(T_c)
    return result
