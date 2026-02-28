from math import *
import numpy as np

def check_harmony(R_0, R_nc, h_c, **kwargs):
    return (R_0) - (R_nc + 1 / h_c)
