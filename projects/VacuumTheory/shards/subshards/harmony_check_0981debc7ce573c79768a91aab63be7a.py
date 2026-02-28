from math import *
import numpy as np

def check_harmony(M, P, T, W, q, **kwargs):
    return (q) - (W * (359 / M) * (760 / P) * (T / 492) * (1/60))
