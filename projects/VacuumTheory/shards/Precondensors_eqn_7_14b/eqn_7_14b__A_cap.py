from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_7_14b__A(self, Q_condensor_heat_duty: float, U: float, del_T_1: float, del_T_2: float, **kwargs):
    # [.pyeqn] A = (Q_condensor_heat_duty / (U * (del_T_1 - del_T_2))) / ln(del_T_1 - del_T_2)
    result = []
    A = Q_condensor_heat_duty/(U*(del_T_1 - del_T_2)*log(del_T_1 - del_T_2))
    result.append(A)
    return result
