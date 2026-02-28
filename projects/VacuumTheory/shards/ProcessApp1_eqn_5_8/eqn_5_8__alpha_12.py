from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_5_8__alpha_12(P_0_1: float, P_0_2: float, gamma_1: float, gamma_2: float, **kwargs):
    # [.pyeqn] alpha_12 = gamma_1 * P_0_1 / (gamma_2 * P_0_2)
    result = []
    alpha_12 = P_0_1*gamma_1/(P_0_2*gamma_2)
    result.append(alpha_12)
    return result
