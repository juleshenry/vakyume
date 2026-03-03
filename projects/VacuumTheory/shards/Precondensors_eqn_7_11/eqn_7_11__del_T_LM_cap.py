from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_7_11__del_T_LM(self, Q_condensor_heat_duty: float, U_v: float, V_c: float, **kwargs):
    # [.pyeqn] V_c = Q_condensor_heat_duty / (U_v * del_T_LM)
    result = []
    del_T_LM = Q_condensor_heat_duty/(U_v*V_c)
    result.append(del_T_LM)
    return result
