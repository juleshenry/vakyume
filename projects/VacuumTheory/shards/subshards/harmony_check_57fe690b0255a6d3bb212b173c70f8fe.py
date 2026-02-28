from math import *
import numpy as np

def check_harmony(M, R, T, g_c, k, v_s, **kwargs):
    return (v_s) - ((k * g_c * R / M * T) ** 0.5)
