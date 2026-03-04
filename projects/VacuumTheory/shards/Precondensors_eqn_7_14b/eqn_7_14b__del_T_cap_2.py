from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_7_14b__del_T_2(
    self, A: float, Q_condensor_heat_duty: float, U: float, del_T_1: float, **kwargs
):
    # [.pyeqn] A = (Q_condensor_heat_duty / (U * (del_T_1 - del_T_2))) / ln(del_T_1 - del_T_2)
    result = []
    del_T_2 = del_T_1 - exp(LambertW(Q_condensor_heat_duty / (A * U)))
    result.append(del_T_2)
    return result
