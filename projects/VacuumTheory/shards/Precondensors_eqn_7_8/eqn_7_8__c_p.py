from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_7_8__c_p(self, L_c: float, Q_condensor_heat_duty: float, del_T: float, **kwargs):
    # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T)
    result = []
    c_p = Q_condensor_heat_duty/(L_c*del_T)
    result.append(c_p)
    return result
