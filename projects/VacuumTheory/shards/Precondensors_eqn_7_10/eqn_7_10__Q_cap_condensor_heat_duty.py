from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_7_10__Q_condensor_heat_duty(self, L_c_P: float, del_T: float, **kwargs):
    # [.pyeqn] L_c_P = Q_condensor_heat_duty / (500 * del_T)
    result = []
    Q_condensor_heat_duty = 500 * L_c_P * del_T
    result.append(Q_condensor_heat_duty)
    return result
