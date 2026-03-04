from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_5_10a__V_1(self, D: float, L_0: float, **kwargs):
    # [.pyeqn] L_0 / V_1 = (L_0 / D) / (L_0 / D + 1)
    # Solved symbolically for V_1
    result = []
    result.append(D + L_0)
    return result
    return [V_1]
