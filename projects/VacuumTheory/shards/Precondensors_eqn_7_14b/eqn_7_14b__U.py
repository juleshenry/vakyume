from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_7_14b__U(
    self,
    A: float,
    Q_condensor_heat_duty: float,
    del_T_1: float,
    del_T_2: float,
    **kwargs,
):
    # [.pyeqn] A = (Q_condensor_heat_duty / (U * (del_T_1 - del_T_2))) / ln(del_T_1 - del_T_2)
    result = []
    U = Q_condensor_heat_duty / (A * (del_T_1 - del_T_2) * log(del_T_1 - del_T_2))
    result.append(U)
    return result
