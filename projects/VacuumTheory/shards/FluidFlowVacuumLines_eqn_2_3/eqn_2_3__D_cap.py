from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_3__D(kn: float, lambd: float, **kwargs):
    # [.pyeqn] kn = lambd / D
    result = []
    D = lambd/kn
    result.append(D)
    return result
