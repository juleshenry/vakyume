from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_10_2__Q_gas(self, PS: float, V: float, dP: float, dt: float, **kwargs):
    # [.pyeqn] PS = - V * dP / dt + Q_gas
    result = []
    Q_gas = PS + V*dP/dt
    result.append(Q_gas)
    return result
