from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_10_2__V(self, PS: float, Q_gas: float, dP: float, dt: float, **kwargs):
    # [.pyeqn] PS = - V * dP / dt + Q_gas
    result = []
    V = dt*(-PS + Q_gas)/dP
    result.append(V)
    return result
