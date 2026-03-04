from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_11_2__SP_1(
    self,
    Q: float,
    Q_0: float,
    Q_external_gas_throughput: float,
    SP_2: float,
    S_vol_pump_speed: float,
    V: float,
    t: float,
    **kwargs,
):
    # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
    result = []
    SP_1 = (
        Q_0
        + Q_external_gas_throughput
        + (-Q - Q_0 + SP_2) * exp(S_vol_pump_speed * t / V)
    )
    result.append(SP_1)
    return result
