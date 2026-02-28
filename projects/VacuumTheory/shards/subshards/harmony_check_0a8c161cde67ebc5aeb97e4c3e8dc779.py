from math import *
import numpy as np

def check_harmony(PS, Q_gas, V, dP, dt, **kwargs):
    return (PS) - (- V * dP / dt + Q_gas)
