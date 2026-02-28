from math import *
import numpy as np

def check_harmony(D, L_0, R, **kwargs):
    return ((L_0 / D) / (L_0 / D + 1)) - (R / (R + 1))
