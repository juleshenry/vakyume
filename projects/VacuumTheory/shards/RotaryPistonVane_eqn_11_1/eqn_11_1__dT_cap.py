from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_11_1__dT(self, PS: float, Q_0: float, Q_external_gas_throughput: float, V: float, dP: float, **kwargs):
    # [.pyeqn] PS = -V * dP / dT + Q_external_gas_throughput + Q_0
    result = []
    dT = V*dP/(-PS + Q_0 + Q_external_gas_throughput)
    result.append(dT)
    return result
