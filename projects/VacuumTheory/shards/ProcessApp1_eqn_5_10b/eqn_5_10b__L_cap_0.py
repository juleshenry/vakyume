from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_5_10b__L_0(self, R: float, V_1: float, **kwargs):
    # [.pyeqn] L_0 / V_1 = R / (R + 1)
    result = []
    L_0 = R * V_1 / (R + 1)
    result.append(L_0)
    return result
