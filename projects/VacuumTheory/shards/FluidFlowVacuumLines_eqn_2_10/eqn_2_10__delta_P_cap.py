from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_2_10__delta_P(self, Suc_Pres: float, oper_press: float, **kwargs):
    # [.pyeqn] Suc_Pres = oper_press - delta_P
    result = []
    delta_P = -Suc_Pres + oper_press
    result.append(delta_P)
    return result
