from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_10_2__Q_gas(self, PS: float, V: float, dP: float, dt: float, **kwargs):
    # [.pyeqn] PS = - V * dP / dt + Q_gas
    result = []
    Q_gas = PS + V*dP/dt
    result.append(Q_gas)
    return result
