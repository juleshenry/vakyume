from math import *
import numpy as np

def check_harmony(A_C, H_2, V_P, **kwargs):
    return (V_P) - (A_C * H_2)
