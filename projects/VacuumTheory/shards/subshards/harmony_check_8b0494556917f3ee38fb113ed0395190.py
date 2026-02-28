from math import *
import numpy as np

def check_harmony(N_mfw, Q_gas, T, **kwargs):
    return (Q_gas) - (9.25 * N_mfw * T)
