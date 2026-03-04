from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_24_9__c(self, K: float, m_0: float, u: float, **kwargs):
    # [.pyeqn] K = (1 - u ** 2 / c ** 2) * m_0 * c ** 2
    result = []
    c = -sqrt((K + m_0 * u**2) / m_0)
    result.append(c)
    c = sqrt((K + m_0 * u**2) / m_0)
    result.append(c)
    return result
