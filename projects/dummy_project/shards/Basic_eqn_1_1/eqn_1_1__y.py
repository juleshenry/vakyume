from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_1_1__y(x: float, z: float, **kwargs):
    # [.pyeqn] z = x + y
    result = []
    y = -x + z
    result.append(y)
    return result
