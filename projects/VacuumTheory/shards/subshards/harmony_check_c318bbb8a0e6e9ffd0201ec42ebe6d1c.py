from math import *
import numpy as np

def check_harmony(P, p_c, p_nc, **kwargs):
    return (p_nc) - (P - p_c)
