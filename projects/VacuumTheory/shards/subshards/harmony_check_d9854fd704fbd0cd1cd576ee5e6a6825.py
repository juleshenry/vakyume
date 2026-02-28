from math import *
import numpy as np

def check_harmony(Q_v, delta_h_v, w_v, **kwargs):
    return (w_v) - (12000 * Q_v / delta_h_v)
