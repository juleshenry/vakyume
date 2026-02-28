from math import *
import numpy as np

def check_harmony(M, P, R, T, V, m, **kwargs):
    return (P * V) - (m / M * R * T)
