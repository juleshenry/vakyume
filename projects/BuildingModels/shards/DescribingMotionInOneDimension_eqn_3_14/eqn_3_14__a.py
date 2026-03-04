from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_3_14__a(self, a_A: float, **kwargs):
    # [.pyeqn] a = a_A
    result = []
    a = a_A
    result.append(a)
    return result
