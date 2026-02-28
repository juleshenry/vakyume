from math import *
import numpy as np

def check_harmony(P, S_0, S_Th, p_0, p_s, **kwargs):
    return (S_Th) - (S_0 * ((P-p_s) / (P - p_0)) ** 0.6)
