from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_10_3__T(self, N_mfw: float, Q_gas: float, **kwargs):
    # [.pyeqn] Q_gas = 9.25 * N_mfw * T
    result = []
    T = 0.108108108108108*Q_gas/N_mfw
    result.append(T)
    return result
