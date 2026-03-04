from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_24_9__m_0(self, K: float, c: float, u: float, **kwargs):
    # [.pyeqn] K = (1 - u ** 2 / c ** 2) * m_0 * c ** 2
    result = []
    m_0 = K / (c**2 - u**2)
    result.append(m_0)
    return result
