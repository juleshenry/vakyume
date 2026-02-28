from math import *
import numpy as np

def check_harmony(H_1, H_2, P, V, V_P, **kwargs):
    return (P) - (V_P * (H_2 - H_1) / (V - V_P))
