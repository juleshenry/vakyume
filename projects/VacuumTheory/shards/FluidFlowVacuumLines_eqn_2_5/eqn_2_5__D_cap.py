from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_5__D(L: float, delta_P: float, mu: float, q: float, **kwargs):
    # [.pyeqn] q = 3.141592653589793 * (D ** 4) * delta_P / (128 * L * mu)
    result = []
    D = -2.52647511098426*I*(L*mu*q/delta_P)**(1/4)
    result.append(D)
    D = 2.52647511098426*I*(L*mu*q/delta_P)**(1/4)
    result.append(D)
    D = -2.52647511098426*(L*mu*q/delta_P)**(1/4)
    result.append(D)
    D = 2.52647511098426*(L*mu*q/delta_P)**(1/4)
    result.append(D)
    return result
