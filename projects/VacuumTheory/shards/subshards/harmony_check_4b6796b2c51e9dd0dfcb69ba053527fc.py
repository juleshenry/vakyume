from math import *
import numpy as np

def check_harmony(H_i, p_i, x_i, **kwargs):
    return (p_i) - (x_i * H_i)
