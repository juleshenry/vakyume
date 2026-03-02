from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_10_3__N_mfw(self, Q_gas: float, T: float, **kwargs):
    # [.pyeqn] Q_gas = 9.25 * N_mfw * T
    result = []
    N_mfw = 0.108108108108108*Q_gas/T
    result.append(N_mfw)
    return result
