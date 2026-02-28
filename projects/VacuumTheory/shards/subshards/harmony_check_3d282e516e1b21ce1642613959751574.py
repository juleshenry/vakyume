from math import *
import numpy as np

def check_harmony(D_eq, R_ll, **kwargs):
    return (D_eq) - (4 * R_ll)
