from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_9_4__r(self, AEL: float, SC: float, w_s: float, **kwargs):
    # [.pyeqn] w_s = AEL * r * SC
    result = []
    r = w_s / (AEL * SC)
    result.append(r)
    return result
