from cmath import *
from math import e, pi
import numpy as np

def check_harmony(T, k, m, v_a, **kwargs):
    return (v_a) - (((8 * k * T) / (3.141592653589793 * m)) ** 0.5)
