from math import *
import numpy as np

def check_harmony(H_1, H_2, P, P_P, **kwargs):
    return (P_P - P) - (H_2 - H_1)
