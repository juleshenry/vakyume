from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_5_8__gamma_1(self, P_0_1: float, P_0_2: float, alpha_12: float, gamma_2: float, **kwargs):
    # [.pyeqn] alpha_12 = gamma_1 * P_0_1 / (gamma_2 * P_0_2)
    result = []
    gamma_1 = P_0_2*alpha_12*gamma_2/P_0_1
    result.append(gamma_1)
    return result
