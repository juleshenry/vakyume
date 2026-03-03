from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_7_14a__U(self, A: float, Q_condensor_heat_duty: float, del_T_LM: float, **kwargs):
    # [.pyeqn] A = Q_condensor_heat_duty / (U * del_T_LM)
    result = []
    U = Q_condensor_heat_duty/(A*del_T_LM)
    result.append(U)
    return result
