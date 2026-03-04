from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_11_4__p_s(self, p_g: float, p_v: float, **kwargs):
    # [.pyeqn] p_v / (p_v + p_g) = p_v / p_s
    result = []
    p_s = p_g + p_v
    result.append(p_s)
    return result
