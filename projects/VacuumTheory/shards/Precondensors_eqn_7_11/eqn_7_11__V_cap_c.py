from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_7_11__V_c(self, Q_condensor_heat_duty: float, U_v: float, del_T_LM: float, **kwargs):
    # [.pyeqn] V_c = Q_condensor_heat_duty / (U_v * del_T_LM)
    result = []
    V_c = Q_condensor_heat_duty/(U_v*del_T_LM)
    result.append(V_c)
    return result
