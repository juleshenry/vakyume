from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_7_9__Q_condensor_heat_duty(
    self, L_c: float, c_p: float, del_T: float, rho: float, **kwargs
):
    # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T * rho * 8.02)
    result = []
    Q_condensor_heat_duty = 8.02 * L_c * c_p * del_T * rho
    result.append(Q_condensor_heat_duty)
    return result
