from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_1_1__y(self, x: float, z: float, **kwargs):
    # [.pyeqn] z = x + y
    result = []
    y = -x + z
    result.append(y)
    return result
