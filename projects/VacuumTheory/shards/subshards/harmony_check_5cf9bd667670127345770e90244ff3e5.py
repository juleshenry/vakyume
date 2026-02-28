from math import *
import numpy as np

def check_harmony(P, S_Th, S_p, p_s, **kwargs):
    return (S_p) - (S_Th * (P - p_s) / P)
