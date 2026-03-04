from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_5_10b__R(self, L_0: float, V_1: float, **kwargs):
    # [.pyeqn] L_0 / V_1 = R / (R + 1)
    result = []
    R = -L_0 / (L_0 - V_1)
    result.append(R)
    return result
