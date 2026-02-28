from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_7_12__U(A: float, Q_condensor_heat_duty: float, del_T: float, **kwargs):
    # [.pyeqn] Q_condensor_heat_duty = U * A * del_T
    result = []
    U = Q_condensor_heat_duty/(A*del_T)
    result.append(U)
    return result
