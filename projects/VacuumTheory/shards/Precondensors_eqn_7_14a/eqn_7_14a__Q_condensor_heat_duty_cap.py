from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_7_14a__Q_condensor_heat_duty(self, A: float, U: float, del_T_LM: float, **kwargs):
    # [.pyeqn] A = Q_condensor_heat_duty / (U * del_T_LM)
    result = []
    Q_condensor_heat_duty = A*U*del_T_LM
    result.append(Q_condensor_heat_duty)
    return result
