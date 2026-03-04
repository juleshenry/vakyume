from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_24_9__K(self, c: float, m_0: float, u: float, **kwargs):
    # [.pyeqn] K = (1 - u ** 2 / c ** 2) * m_0 * c ** 2
    result = []
    K = m_0 * (c**2 - u**2)
    result.append(K)
    return result
