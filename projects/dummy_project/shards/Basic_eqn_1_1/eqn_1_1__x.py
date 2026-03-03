from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_1_1__x(self, y: float, z: float, **kwargs):
    # [.pyeqn] z = x + y
    result = []
    x = -y + z
    result.append(x)
    return result
