from math import *
import numpy as np

def check_harmony(p_g, p_s, p_v, **kwargs):
    return (p_v / (p_v + p_g)) - (p_v / p_s)
