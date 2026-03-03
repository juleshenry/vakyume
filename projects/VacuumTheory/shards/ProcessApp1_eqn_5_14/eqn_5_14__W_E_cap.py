from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_5_14__W_E(self, M: float, P_0: float, T: float, **kwargs):
    # [.pyeqn] W_E = 0.0583 * P_0 * (M / T) ** 0.5
    result = []
    W_E = 0.0583*P_0*sqrt(M/T)
    result.append(W_E)
    return result
