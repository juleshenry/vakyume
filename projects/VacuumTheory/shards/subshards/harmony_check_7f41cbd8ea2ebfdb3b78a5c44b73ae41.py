from math import *
import numpy as np

def check_harmony(P, P_P, V, V_P, **kwargs):
    return (P_P) - (P * (V / V_P))
