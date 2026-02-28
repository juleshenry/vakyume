from math import *
import numpy as np

def check_harmony(PS, Q_0, Q_external_gas_throughput, V, dP, dT, **kwargs):
    return (PS) - (-V * dP / dT + Q_external_gas_throughput + Q_0)
