from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_11_3__t(self, F_s: float, t_c: float, **kwargs):
    # [.pyeqn] t = t_c * F_s
    result = []
    t = F_s * t_c
    result.append(t)
    return result
