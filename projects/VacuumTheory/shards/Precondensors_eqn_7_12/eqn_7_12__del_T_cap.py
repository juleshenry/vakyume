from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_7_12__del_T(A: float, Q_condensor_heat_duty: float, U: float, **kwargs):
    # [.pyeqn] Q_condensor_heat_duty = U * A * del_T
    result = []
    del_T = Q_condensor_heat_duty/(A*U)
    result.append(del_T)
    return result
