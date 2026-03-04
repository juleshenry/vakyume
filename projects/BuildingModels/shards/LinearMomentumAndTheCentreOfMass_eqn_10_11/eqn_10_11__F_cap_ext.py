from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_10_11__F_ext(self, M: float, a_CM: float, **kwargs):
    # [.pyeqn] F_ext = M * a_CM
    result = []
    F_ext = M * a_CM
    result.append(F_ext)
    return result
