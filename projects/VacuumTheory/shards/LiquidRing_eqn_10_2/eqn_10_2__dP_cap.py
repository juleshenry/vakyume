from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_10_2__dP(PS: float, Q_gas: float, V: float, dt: float, **kwargs):
    # [.pyeqn] PS = - V * dP / dt + Q_gas
    result = []
    dP = dt*(-PS + Q_gas)/V
    result.append(dP)
    return result
