from math import *
import numpy as np

def check_harmony(T, V, del_P, leakage, t, **kwargs):
    return (leakage) - (0.0059 * V * del_P / t * 530 / T)
