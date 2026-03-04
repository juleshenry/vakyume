from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_5_10c__L_0(self, D: float, R: float, **kwargs):
    # [.pyeqn] (L_0 / D) / (L_0 / D + 1) = R / (R + 1)
    result = []
    L_0 = D * R
    result.append(L_0)
    return result
