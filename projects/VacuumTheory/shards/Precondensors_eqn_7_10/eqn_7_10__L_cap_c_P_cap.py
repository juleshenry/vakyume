from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_7_10__L_c_P(self, Q_condensor_heat_duty: float, del_T: float, **kwargs):
    # [.pyeqn] L_c_P = Q_condensor_heat_duty / (500 * del_T)
    result = []
    L_c_P = Q_condensor_heat_duty/(500*del_T)
    result.append(L_c_P)
    return result
