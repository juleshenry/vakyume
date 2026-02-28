from math import *
import numpy as np

def check_harmony(P_0_i, p_i, x_i, **kwargs):
    return (p_i) - (x_i * P_0_i)
