from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_7_4ab__P_c(p: float, p_i: float, p_nc: float, **kwargs):
    # [.pyeqn] p_i / p_nc = p_i / (p - P_c)
    result = []
    P_c = p - p_nc
    result.append(P_c)
    return result
