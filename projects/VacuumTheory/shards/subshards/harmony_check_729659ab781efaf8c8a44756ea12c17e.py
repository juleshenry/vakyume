from cmath import *
from math import e, pi
import numpy as np

def check_harmony(P_0_i, gamma_i, p_i, x_i, **kwargs):
    return (p_i) - (x_i * gamma_i * P_0_i)
