from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_3__lambd(self, D: float, kn: float, **kwargs):
    # [.pyeqn] kn = lambd / D
    result = []
    lambd = D*kn
    result.append(lambd)
    return result
