from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_5_13__HETP(self, H_p: float, N_ES: float, **kwargs):
    # [.pyeqn] H_p = N_ES * HETP
    result = []
    HETP = H_p/N_ES
    result.append(HETP)
    return result
