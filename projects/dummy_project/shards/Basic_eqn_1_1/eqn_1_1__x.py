from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_1_1__x(y: float, z: float, **kwargs):
    # [.pyeqn] z = x + y
    result = []
    x = -y + z
    result.append(x)
    return result
