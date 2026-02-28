from math import *
import numpy as np

def check_harmony(T_c, T_s, delta_T, **kwargs):
    return (T_c) - (T_s + delta_T)
