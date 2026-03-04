from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_7_4ab__p_nc(self, P_c: float, p: float, p_i: float, **kwargs):
    # [.pyeqn] p_i / p_nc = p_i / (p - P_c)
    result = []
    p_nc = -P_c + p
    result.append(p_nc)
    return result
