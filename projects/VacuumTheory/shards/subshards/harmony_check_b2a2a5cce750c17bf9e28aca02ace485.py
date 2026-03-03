from cmath import *
from math import e, pi
import numpy as np

def check_harmony(A_d, delta_T, delta_h_i, delta_m, h_d, m_b, t_R, **kwargs):
    return (t_R) - (delta_h_i * m_b * delta_m / (A_d * h_d * delta_T))
