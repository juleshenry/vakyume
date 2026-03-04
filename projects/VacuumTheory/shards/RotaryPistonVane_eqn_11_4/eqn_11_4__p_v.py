from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_11_4__p_v(self, p_g: float, p_s: float, **kwargs):
    # [.pyeqn] p_v / (p_v + p_g) = p_v / p_s
    result = []
    p_v = 0
    result.append(p_v)
    p_v = -p_g + p_s
    result.append(p_v)
    return result
