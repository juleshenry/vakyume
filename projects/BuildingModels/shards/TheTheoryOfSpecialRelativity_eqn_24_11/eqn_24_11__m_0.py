from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_24_11__m_0(self, E: float, c: float, p: float, **kwargs):
    # [.pyeqn] E ** 2 = p ** 2 * c ** 2 + m_0 ** 2
    result = []
    m_0 = -sqrt((E - c * p) * (E + c * p))
    result.append(m_0)
    m_0 = sqrt((E - c * p) * (E + c * p))
    result.append(m_0)
    return result
