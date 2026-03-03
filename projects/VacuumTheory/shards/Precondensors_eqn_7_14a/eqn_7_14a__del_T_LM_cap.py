from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_7_14a__del_T_LM(self, A: float, Q_condensor_heat_duty: float, U: float, **kwargs):
    # [.pyeqn] A = Q_condensor_heat_duty / (U * del_T_LM)
    result = []
    del_T_LM = Q_condensor_heat_duty/(A*U)
    result.append(del_T_LM)
    return result
