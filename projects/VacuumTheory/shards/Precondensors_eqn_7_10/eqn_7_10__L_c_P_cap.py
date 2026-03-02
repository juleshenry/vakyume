from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_7_10__L_c_P(self, Q_condensor_heat_duty: float, del_T: float, **kwargs):
    # [.pyeqn] L_c_P = Q_condensor_heat_duty / (500 * del_T)
    result = []
    L_c_P = Q_condensor_heat_duty/(500*del_T)
    result.append(L_c_P)
    return result
