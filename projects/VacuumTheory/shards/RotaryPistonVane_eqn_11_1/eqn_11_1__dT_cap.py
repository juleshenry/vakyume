from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_11_1__dT(self, PS: float, Q_0: float, Q_external_gas_throughput: float, V: float, dP: float, **kwargs):
    # [.pyeqn] PS = -V * dP / dT + Q_external_gas_throughput + Q_0
    result = []
    dT = V*dP/(-PS + Q_0 + Q_external_gas_throughput)
    result.append(dT)
    return result
