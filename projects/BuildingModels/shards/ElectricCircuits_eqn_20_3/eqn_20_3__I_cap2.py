from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_20_3__I2(self, I1: float, I3: float, **kwargs):
    # [.pyeqn] I1 = I2 + I3
    result = []
    I2 = I1 - I3
    result.append(I2)
    return result
