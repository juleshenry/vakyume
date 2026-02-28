from math import *
import numpy as np

def check_harmony(M, P, P_i_0, W_air, W_i, epsilon_i, p_c, x_i, **kwargs):
    return (W_i) - (W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c)))
