from math import *
import numpy as np

def check_harmony(C_L, C_T, F_p, **kwargs):
    return (C_T) - (C_L * F_p)
