from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_34__C_1(self, C: float, C_2: float, D: float, L: float, P_p: float, mu: float, **kwargs):
    # [.pyeqn] C = C_1 * (D ** 4 / (mu * L)) * P_p + C_2 * (D ** 3 / L)
    result = []
    C_1 = mu*(C*L - C_2*D**3)/(D**4*P_p)
    result.append(C_1)
    return result
