from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_3_4__V(KAPPA: float, P: float, **kwargs):
    # [.pyeqn] P * V = KAPPA
    result = []
    V = KAPPA/P
    result.append(V)
    return result
