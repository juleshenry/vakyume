from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_11_3__t_c(self, F_s: float, t: float, **kwargs):
    # [.pyeqn] t = t_c * F_s
    result = []
    t_c = t / F_s
    result.append(t_c)
    return result
