from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_11_2__S_vol_pump_speed(Q: float, Q_0: float, Q_external_gas_throughput: float, SP_1: float, SP_2: float, V: float, t: float, **kwargs):
    # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
    result = []
    S_vol_pump_speed = V*log((Q_0 + Q_external_gas_throughput - SP_1)/(Q + Q_0 - SP_2))/t
    result.append(S_vol_pump_speed)
    return result
