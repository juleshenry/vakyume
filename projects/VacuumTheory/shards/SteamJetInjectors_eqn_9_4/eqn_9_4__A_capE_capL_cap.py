from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_9_4__AEL(self, SC: float, r: float, w_s: float, **kwargs):
    # [.pyeqn] w_s = AEL * r * SC
    result = []
    AEL = w_s / (SC * r)
    result.append(AEL)
    return result
