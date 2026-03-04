from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_11_2__SP_2(self, Q: float, Q_0: float, Q_external_gas_throughput: float, SP_1: float, S_vol_pump_speed: float, V: float, t: float, **kwargs):
    # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
    result = []
    SP_2 = (-Q_0 - Q_external_gas_throughput + SP_1 + (Q + Q_0)*exp(S_vol_pump_speed*t/V))*exp(-S_vol_pump_speed*t/V)
    result.append(SP_2)
    return result
