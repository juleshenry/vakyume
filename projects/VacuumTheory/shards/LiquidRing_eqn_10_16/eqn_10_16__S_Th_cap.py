from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_10_16__S_Th(self, P: float, S_0: float, p_0: float, **kwargs):
    # [.pyeqn] S_Th = S_0 * (P / (P - p_0)) ** 0.6
    result = []
    S_Th = S_0 * (P / (P - p_0)) ** (3 / 5)
    result.append(S_Th)
    return result
