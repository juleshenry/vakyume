from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_11_1__V(self, PS: float, Q_0: float, Q_external_gas_throughput: float, dP: float, dT: float, **kwargs):
    # [.pyeqn] PS = -V * dP / dT + Q_external_gas_throughput + Q_0
    result = []
    V = dT*(-PS + Q_0 + Q_external_gas_throughput)/dP
    result.append(V)
    return result
