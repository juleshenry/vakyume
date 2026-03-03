from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_7_4ac__p_i(self, P_c: float, n_i: float, n_nc: float, p: float, **kwargs):
    # [.pyeqn] n_i / n_nc = p_i / (p - P_c)
    result = []
    p_i = n_i*(-P_c + p)/n_nc
    result.append(p_i)
    return result
