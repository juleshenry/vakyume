from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_7_4a__p_c(self, P: float, p_nc: float, **kwargs):
    # [.pyeqn] p_nc = P - p_c
    result = []
    p_c = P - p_nc
    result.append(p_c)
    return result
