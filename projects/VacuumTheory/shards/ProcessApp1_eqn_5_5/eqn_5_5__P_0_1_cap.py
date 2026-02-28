from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_5_5__P_0_1(P_0_2: float, alpha_12: float, **kwargs):
    # [.pyeqn] alpha_12 = P_0_1 / P_0_2
    result = []
    P_0_1 = P_0_2*alpha_12
    result.append(P_0_1)
    return result
