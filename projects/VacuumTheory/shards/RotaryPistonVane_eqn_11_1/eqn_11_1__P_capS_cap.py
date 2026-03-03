from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_11_1__PS(self, Q_0: float, Q_external_gas_throughput: float, V: float, dP: float, dT: float, **kwargs):
    # [.pyeqn] PS = -V * dP / dT + Q_external_gas_throughput + Q_0
    result = []
    PS = Q_0 + Q_external_gas_throughput - V*dP/dT
    result.append(PS)
    return result
