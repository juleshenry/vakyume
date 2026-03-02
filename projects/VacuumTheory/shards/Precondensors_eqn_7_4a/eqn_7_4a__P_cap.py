from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_7_4a__P(self, p_c: float, p_nc: float, **kwargs):
    # [.pyeqn] p_nc = P - p_c
    result = []
    P = p_c + p_nc
    result.append(P)
    return result
