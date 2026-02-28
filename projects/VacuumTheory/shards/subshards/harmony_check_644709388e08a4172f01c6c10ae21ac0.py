from math import *
import numpy as np

def check_harmony(K_i, x_i, y_i, **kwargs):
    return (K_i) - (y_i / x_i)
