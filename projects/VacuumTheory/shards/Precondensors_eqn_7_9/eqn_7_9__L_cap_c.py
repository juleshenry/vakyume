from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_7_9__L_c(
    self, Q_condensor_heat_duty: float, c_p: float, del_T: float, rho: float, **kwargs
):
    # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T * rho * 8.02)
    result = []
    L_c = 0.124688279301746 * Q_condensor_heat_duty / (c_p * del_T * rho)
    result.append(L_c)
    return result
