from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_7_8__del_T(self, L_c: float, Q_condensor_heat_duty: float, c_p: float, **kwargs):
    # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T)
    result = []
    del_T = Q_condensor_heat_duty/(L_c*c_p)
    result.append(del_T)
    return result
