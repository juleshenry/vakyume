from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_2_8__T_c(self, M: float, P_c: float, mu_c: float, **kwargs):
    # [.pyeqn] mu_c = (7.7 * (M ** 0.5) * P_c ** (2 / 3)) / T_c ** (1 / 6)
    result = []
    T_c = 208422.380089 * M**3 * P_c**4 / mu_c**6
    result.append(T_c)
    return result
