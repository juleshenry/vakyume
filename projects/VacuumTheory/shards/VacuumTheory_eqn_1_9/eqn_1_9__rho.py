from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_1_9__rho(M: float, P: float, R: float, T: float, **kwargs):
    # [.pyeqn] rho = P * M / (R * T)
    result = []
    rho = M*P/(R*T)
    result.append(rho)
    return result
