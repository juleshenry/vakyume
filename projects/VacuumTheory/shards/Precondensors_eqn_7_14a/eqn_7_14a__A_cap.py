from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_7_14a__A(self, Q_condensor_heat_duty: float, U: float, del_T_LM: float, **kwargs):
    # [.pyeqn] A = Q_condensor_heat_duty / (U * del_T_LM)
    result = []
    A = Q_condensor_heat_duty/(U*del_T_LM)
    result.append(A)
    return result
