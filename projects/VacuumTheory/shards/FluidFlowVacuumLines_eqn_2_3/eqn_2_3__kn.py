from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_3__kn(self, D: float, lambd: float, **kwargs):
    # [.pyeqn] kn = lambd / D
    result = []
    kn = lambd/D
    result.append(kn)
    return result
