from math import *
import numpy as np

def check_harmony(R, T, V, n, p, **kwargs):
    return (p * V) - (n * R * T)
