from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_20_3__I3(self, I1: float, I2: float, **kwargs):
    # [.pyeqn] I1 = I2 + I3
    result = []
    I3 = I1 - I2
    result.append(I3)
    return result
