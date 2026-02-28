from math import *
import numpy as np

def check_harmony(D_0, D_LM, D_i, R_f_0, R_fi, U_0, h_0, h_i, k_w, x_w, **kwargs):
    return (1 / U_0) - (1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i))
