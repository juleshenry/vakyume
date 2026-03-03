from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_10_2__dt(self, PS: float, Q_gas: float, V: float, dP: float, **kwargs):
    # [.pyeqn] PS = - V * dP / dt + Q_gas
    result = []
    dt = -V * dP / (PS - Q_gas)
    result.append(dt)
    return result
