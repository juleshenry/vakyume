from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_5_14__M(self, P_0: float, T: float, W_E: float, **kwargs):
    # [.pyeqn] W_E = 0.0583 * P_0 * (M / T) ** 0.5
    result = []
    M = 294.213699178261*T*W_E**2/P_0**2
    result.append(M)
    return result
