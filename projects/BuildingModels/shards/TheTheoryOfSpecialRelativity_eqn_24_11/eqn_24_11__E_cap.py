from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_24_11__E(self, c: float, m_0: float, p: float, **kwargs):
    # [.pyeqn] E ** 2 = p ** 2 * c ** 2 + m_0 ** 2
    result = []
    E = -sqrt(c**2 * p**2 + m_0**2)
    result.append(E)
    E = sqrt(c**2 * p**2 + m_0**2)
    result.append(E)
    return result
