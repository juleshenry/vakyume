from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_2_10__oper_press(self, Suc_Pres: float, delta_P: float, **kwargs):
    # [.pyeqn] Suc_Pres = oper_press - delta_P
    result = []
    oper_press = Suc_Pres + delta_P
    result.append(oper_press)
    return result
