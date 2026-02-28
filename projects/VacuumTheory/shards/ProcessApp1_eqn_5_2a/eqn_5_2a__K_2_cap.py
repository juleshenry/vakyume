from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_5_2a__K_2(K_1: float, alpha_1_2: float, **kwargs):
    # [.pyeqn] alpha_1_2 = K_1 / K_2
    result = []
    K_2 = K_1/alpha_1_2
    result.append(K_2)
    return result
