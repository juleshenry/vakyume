from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_7_8__c_p(L_c: float, Q_condensor_heat_duty: float, del_T: float, **kwargs):
    # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T)
    result = []
    c_p = Q_condensor_heat_duty/(L_c*del_T)
    result.append(c_p)
    return result
