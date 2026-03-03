from cmath import *
from math import e, pi
import numpy as np

def check_harmony(N_i, N_nc, P, P_c, p_i, **kwargs):
    return (N_i) - (N_nc * (p_i) / (P - P_c))
