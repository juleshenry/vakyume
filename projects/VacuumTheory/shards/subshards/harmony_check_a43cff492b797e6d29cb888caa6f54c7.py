from math import *
import numpy as np

def check_harmony(AEL, SC, r, w_s, **kwargs):
    return (w_s) - (AEL * r * SC)
