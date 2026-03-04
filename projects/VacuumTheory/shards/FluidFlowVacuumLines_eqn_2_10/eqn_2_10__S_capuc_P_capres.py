from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_2_10__Suc_Pres(self, delta_P: float, oper_press: float, **kwargs):
    # [.pyeqn] Suc_Pres = oper_press - delta_P
    result = []
    Suc_Pres = -delta_P + oper_press
    result.append(Suc_Pres)
    return result
