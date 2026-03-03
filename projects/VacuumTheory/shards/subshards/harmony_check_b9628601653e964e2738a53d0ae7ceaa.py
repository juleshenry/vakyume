from cmath import *
from math import e, pi
import numpy as np

def check_harmony(A, C, F_t, M, T, **kwargs):
    return (C) - (38.3 * (T * A * F_t / M) ** 0.5)
