from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_7_12__Q_condensor_heat_duty(self, A: float, U: float, del_T: float, **kwargs):
    # [.pyeqn] Q_condensor_heat_duty = U * A * del_T
    result = []
    Q_condensor_heat_duty = A*U*del_T
    result.append(Q_condensor_heat_duty)
    return result
