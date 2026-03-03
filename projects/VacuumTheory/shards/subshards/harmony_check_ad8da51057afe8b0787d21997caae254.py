from cmath import *
from math import e, pi
import numpy as np

def check_harmony(M, P_0, T, W_E, **kwargs):
    return (W_E) - (0.0583 * P_0 * (M / T) ** 0.5)
