from math import *
import numpy as np

def check_harmony(W, W_T, sum_individual_leak_rates, **kwargs):
    return (W_T) - (W + sum_individual_leak_rates)
