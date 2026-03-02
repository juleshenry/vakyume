from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_7_4ab__p_i(self, P_c: float, p: float, p_nc: float, **kwargs):
    # [.pyeqn] p_i / p_nc = p_i / (p - P_c)
    result = []
    p_i = 0
    result.append(p_i)
    return result
