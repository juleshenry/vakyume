from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_2_8__P_c(self, M: float, T_c: float, mu_c: float, **kwargs):
    # [.pyeqn] mu_c = (7.7 * (M ** 0.5) * P_c ** (2 / 3)) / T_c ** (1 / 6)
    result = []
    P_c = -0.046801946114055 * (T_c**0.166666666666667 * mu_c / M**0.5) ** (3 / 2)
    result.append(P_c)
    P_c = 0.046801946114055 * (T_c**0.166666666666667 * mu_c / M**0.5) ** (3 / 2)
    result.append(P_c)
    return result
