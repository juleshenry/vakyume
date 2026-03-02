from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_11_1__Q_0(self, PS: float, Q_external_gas_throughput: float, V: float, dP: float, dT: float, **kwargs):
    # [.pyeqn] PS = -V * dP / dT + Q_external_gas_throughput + Q_0
    result = []
    Q_0 = PS - Q_external_gas_throughput + V*dP/dT
    result.append(Q_0)
    return result
