from math import *
import numpy as np

def check_harmony(P_m, d_n, rho_s, w_s, **kwargs):
    return (w_s) - (865.8 * d_n**2 * (P_m * rho_s) ** 0.5)
