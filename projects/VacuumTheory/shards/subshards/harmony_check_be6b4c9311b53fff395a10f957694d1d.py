from math import *
import numpy as np

def check_harmony(n_i, n_nc, p_i, p_nc, **kwargs):
    return (n_i / n_nc) - (p_i / p_nc)
