from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_7_11__U_v(self, Q_condensor_heat_duty: float, V_c: float, del_T_LM: float, **kwargs):
    # [.pyeqn] V_c = Q_condensor_heat_duty / (U_v * del_T_LM)
    result = []
    U_v = Q_condensor_heat_duty/(V_c*del_T_LM)
    result.append(U_v)
    return result
