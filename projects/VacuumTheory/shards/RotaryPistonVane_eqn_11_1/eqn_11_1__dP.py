from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_11_1__dP(self, PS: float, Q_0: float, Q_external_gas_throughput: float, V: float, dT: float, **kwargs):
    # [.pyeqn] PS = -V * dP / dT + Q_external_gas_throughput + Q_0
    result = []
    dP = dT*(-PS + Q_0 + Q_external_gas_throughput)/V
    result.append(dP)
    return result
