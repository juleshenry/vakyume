from cmath import *
from math import e, pi
import numpy as np

def check_harmony(P_i_0, p_i, x_i, **kwargs):
    return (p_i) - (x_i * P_i_0)
