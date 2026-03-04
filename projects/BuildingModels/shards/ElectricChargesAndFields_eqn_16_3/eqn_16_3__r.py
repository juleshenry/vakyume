from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_16_3__r(self, F_g: float, G: float, m_e: float, m_p: float, **kwargs):
    # [.pyeqn] F_g = G * m_e * m_p / r ** 2
    result = []
    r = -sqrt(G * m_e * m_p / F_g)
    result.append(r)
    r = sqrt(G * m_e * m_p / F_g)
    result.append(r)
    return result
