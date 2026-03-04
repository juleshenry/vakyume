from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_9_3__w_j(self, P_s: float, V: float, t_e: float, **kwargs):
    # [.pyeqn] t_e = (2.3 - 0.003 * P_s) * V / w_j
    result = []
    w_j = 0.001*V*(2300.0 - 3.0*P_s)/t_e
    result.append(w_j)
    return result
