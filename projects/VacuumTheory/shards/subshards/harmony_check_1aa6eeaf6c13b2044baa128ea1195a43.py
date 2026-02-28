from math import *
import numpy as np

def check_harmony(L, sum_equivalent_length, sum_pipe, **kwargs):
    return (L) - (sum_pipe + sum_equivalent_length)
