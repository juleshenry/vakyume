from math import *
import numpy as np

def check_harmony(D, L, f, g_c, h_r, v, **kwargs):
    return (h_r) - (f * L * v ** 2 / (D * 2 * g_c))
