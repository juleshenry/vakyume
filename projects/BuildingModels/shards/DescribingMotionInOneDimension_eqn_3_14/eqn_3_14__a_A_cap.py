from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_3_14__a_A(self, a: float, **kwargs):
    # [.pyeqn] a = a_A
    result = []
    a_A = a
    result.append(a_A)
    return result
