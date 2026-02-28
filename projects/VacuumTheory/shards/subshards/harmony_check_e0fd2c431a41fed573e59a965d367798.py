from math import *
import numpy as np

def check_harmony(P, P_0_i, gamma_i, x_i, y_i, **kwargs):
    return (y_i * P) - (x_i * gamma_i * P_0_i)
