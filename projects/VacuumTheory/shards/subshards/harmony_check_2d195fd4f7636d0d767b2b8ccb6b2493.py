from math import *
import numpy as np

def check_harmony(H_2_1, H_2_3, H_2_mi, x_1, x_3, **kwargs):
    return (log(H_2_mi)) - (x_1 * log(H_2_1) + x_3 * log(H_2_3))
