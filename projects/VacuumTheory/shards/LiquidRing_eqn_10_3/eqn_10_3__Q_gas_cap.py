from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_10_3__Q_gas(N_mfw: float, T: float, **kwargs):
    # [.pyeqn] Q_gas = 9.25 * N_mfw * T
    result = []
    Q_gas = 9.25*N_mfw*T
    result.append(Q_gas)
    return result
