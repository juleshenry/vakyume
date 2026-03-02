from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_10__delta_P(self, Suc_Pres: float, oper_press: float, **kwargs):
    # [.pyeqn] Suc_Pres = oper_press - delta_P
    result = []
    delta_P = -Suc_Pres + oper_press
    result.append(delta_P)
    return result
