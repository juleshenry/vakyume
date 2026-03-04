from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_24_6__c(self, u: float, ux: float, v: float, **kwargs):
    # [.pyeqn] ux = u + v / (1 + (v * u) / c ** 2)
    result = []
    c = sqrt(u * v * (-u + ux) / (u - ux + v))
    result.append(c)
    c = -sqrt(-u * v * (u - ux) / (u - ux + v))
    result.append(c)
    return result
