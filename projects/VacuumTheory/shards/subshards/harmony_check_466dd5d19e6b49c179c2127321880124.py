from math import *
import numpy as np

def check_harmony(C_series, geometric_sum_C, **kwargs):
    return (1 / C_series) - (geometric_sum_C)
