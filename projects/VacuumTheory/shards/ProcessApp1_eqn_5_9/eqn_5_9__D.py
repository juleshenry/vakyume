from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_5_9__D(self, L_0, V_1, **kwargs):
    # [.pyeqn] L_0 / V_1 = L_0 / (L_0 + D)
    # Solved symbolically for D
    result = []
    result.append(-L_0 + V_1)
    return result
