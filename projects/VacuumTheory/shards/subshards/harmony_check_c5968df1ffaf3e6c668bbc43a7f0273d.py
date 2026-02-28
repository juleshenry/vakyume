from math import *
import numpy as np

def check_harmony(HETP, H_p, N_ES, **kwargs):
    return (H_p) - (N_ES * HETP)
