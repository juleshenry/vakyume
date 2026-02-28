from math import *
import numpy as np

def check_harmony(A_C, H_2, P, V, **kwargs):
    return (P) - (A_C / V * (H_2) ** 2)
