from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_7_8__L_c(Q_condensor_heat_duty: float, c_p: float, del_T: float, **kwargs):
    # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T)
    result = []
    L_c = Q_condensor_heat_duty/(c_p*del_T)
    result.append(L_c)
    return result
