from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_24_11__c(self, E: float, m_0: float, p: float, **kwargs):
    # [.pyeqn] E ** 2 = p ** 2 * c ** 2 + m_0 ** 2
    result = []
    c = -sqrt((E - m_0) * (E + m_0)) / p
    result.append(c)
    c = sqrt((E - m_0) * (E + m_0)) / p
    result.append(c)
    return result
