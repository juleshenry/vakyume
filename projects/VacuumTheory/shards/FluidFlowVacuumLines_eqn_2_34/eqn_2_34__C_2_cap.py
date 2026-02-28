from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_34__C_2(C: float, C_1: float, D: float, L: float, P_p: float, mu: float, **kwargs):
    # [.pyeqn] C = C_1 * (D ** 4 / (mu * L)) * P_p + C_2 * (D ** 3 / L)
    result = []
    C_2 = C*L/D**3 - C_1*D*P_p/mu
    result.append(C_2)
    return result
