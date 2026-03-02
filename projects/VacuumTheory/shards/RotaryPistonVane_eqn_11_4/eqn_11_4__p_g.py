from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_11_4__p_g(self, p_s: float, p_v: float, **kwargs):
    # [.pyeqn] p_v / (p_v + p_g) = p_v / p_s
    result = []
    p_g = p_s - p_v
    result.append(p_g)
    return result
