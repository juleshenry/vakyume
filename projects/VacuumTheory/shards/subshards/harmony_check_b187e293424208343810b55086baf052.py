from cmath import *
from math import e, pi
import numpy as np


def check_harmony(
    Q, Q_0, Q_external_gas_throughput, SP_1, SP_2, S_vol_pump_speed, V, t, **kwargs
):
    return (t) - (
        V
        / S_vol_pump_speed
        * log((SP_1 - (Q_external_gas_throughput + Q_0)) / (SP_2 - (Q + Q_0)))
    )
