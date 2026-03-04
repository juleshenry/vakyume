from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_5_10a__D(self, L_0: float, V_1: float, **kwargs):
    # [.pyeqn] L_0 / V_1 = (L_0 / D) / (L_0 / D + 1)
    # Solved symbolically for D
    result = []
    result.append(-L_0 + V_1)
    return result
