from math import *
import numpy as np

def check_harmony(P_1, P_2, S_a, V, t, **kwargs):
    return (S_a) - (V / t * log(P_1 / P_2))
