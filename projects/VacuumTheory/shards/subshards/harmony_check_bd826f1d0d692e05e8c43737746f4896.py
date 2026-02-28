from math import *
import numpy as np

def check_harmony(D, L_0, V_1, **kwargs):
    return (L_0 / V_1) - (L_0 / (L_0 + D))
