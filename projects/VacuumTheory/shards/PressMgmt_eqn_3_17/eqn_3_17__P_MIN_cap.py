from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_3_17__P_MIN(self, **kwargs):
    # [.pyeqn] P_MIN = (3.141592653589793 / 4) / (200000)
    result = []
    P_MIN = 0.00000392699081698724
    result.append(P_MIN)
    return result
