from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_7_12__A(self, Q_condensor_heat_duty: float, U: float, del_T: float, **kwargs):
    # [.pyeqn] Q_condensor_heat_duty = U * A * del_T
    result = []
    A = Q_condensor_heat_duty / (U * del_T)
    result.append(A)
    return result
