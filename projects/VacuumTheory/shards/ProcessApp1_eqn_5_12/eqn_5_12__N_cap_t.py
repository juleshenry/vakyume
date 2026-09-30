from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_5_12__N_t(self, Eff: float, N_ES: float, T: float, **kwargs):
    # [.pyeqn] N_t = N_ES / Eff ** T
    result = []
    N_t = N_ES / Eff**T
    result.append(N_t)
    return result
