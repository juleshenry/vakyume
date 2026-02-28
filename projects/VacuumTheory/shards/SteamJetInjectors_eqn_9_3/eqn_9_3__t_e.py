from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_9_3__t_e(P_s: float, V: float, w_j: float, **kwargs):
    # [.pyeqn] t_e = (2.3 - 0.003 * P_s) * V / w_j
    result = []
    t_e = 0.001*V*(2300.0 - 3.0*P_s)/w_j
    result.append(t_e)
    return result
