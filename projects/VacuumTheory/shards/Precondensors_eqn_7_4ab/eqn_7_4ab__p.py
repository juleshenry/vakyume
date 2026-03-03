from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_7_4ab__p(self, P_c: float, p_i: float, p_nc: float, **kwargs):
    # [.pyeqn] p_i / p_nc = p_i / (p - P_c)
    result = []
    p = P_c + p_nc
    result.append(p)
    return result
