from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_7_11__Q_condensor_heat_duty(self, U_v: float, V_c: float, del_T_LM: float, **kwargs):
    # [.pyeqn] V_c = Q_condensor_heat_duty / (U_v * del_T_LM)
    result = []
    Q_condensor_heat_duty = U_v*V_c*del_T_LM
    result.append(Q_condensor_heat_duty)
    return result
