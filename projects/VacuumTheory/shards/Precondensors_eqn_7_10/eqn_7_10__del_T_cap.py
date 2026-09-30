from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_7_10__del_T(self, L_c_P: float, Q_condensor_heat_duty: float, **kwargs):
    # [.pyeqn] L_c_P = Q_condensor_heat_duty / (500 * del_T)
    result = []
    del_T = Q_condensor_heat_duty / (500 * L_c_P)
    result.append(del_T)
    return result
