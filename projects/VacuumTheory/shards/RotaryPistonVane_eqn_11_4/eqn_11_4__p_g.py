from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_11_4__p_g(self, p_s: float, p_v: float, **kwargs):
    # [.pyeqn] p_v / (p_v + p_g) = p_v / p_s
    result = []
    p_g = p_s - p_v
    result.append(p_g)
    return result
