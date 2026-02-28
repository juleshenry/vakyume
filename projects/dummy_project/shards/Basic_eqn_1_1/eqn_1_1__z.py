from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_1_1__z(x: float, y: float, **kwargs):
    # [.pyeqn] z = x + y
    result = []
    z = x + y
    result.append(z)
    return result
