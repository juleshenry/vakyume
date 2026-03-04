from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_10_3__Q_gas(self, N_mfw: float, T: float, **kwargs):
    # [.pyeqn] Q_gas = 9.25 * N_mfw * T
    result = []
    Q_gas = 9.25*N_mfw*T
    result.append(Q_gas)
    return result
