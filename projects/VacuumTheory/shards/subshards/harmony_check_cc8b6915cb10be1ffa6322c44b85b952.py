from cmath import *
from math import e, pi
import numpy as np


def check_harmony(D_0, D_LM, D_i, R_fi, R_fo, R_nc, U_0, h_c, h_i, k_w, x_w, **kwargs):
    return (1 / U_0) - (
        R_nc
        + 1 / h_c
        + R_fo
        + (x_w * D_0) / (k_w * D_LM)
        + R_fi * D_0 / D_i
        + D_0 / (h_i * D_i)
    )
