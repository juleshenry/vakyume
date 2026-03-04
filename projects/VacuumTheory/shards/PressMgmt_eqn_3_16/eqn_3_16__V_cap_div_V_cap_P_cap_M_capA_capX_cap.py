from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_3_16__V_div_V_P_MAX(self, **kwargs):
    # [.pyeqn] V_div_V_P_MAX = 200000 / (3.141592653589793 / 4)
    result = []
    V_div_V_P_MAX = 254647.908947033
    result.append(V_div_V_P_MAX)
    return result
