from cmath import *
from math import e, pi
import numpy as np

def check_harmony(C, C_0, F_t, **kwargs):
    return (C) - (C_0 * F_t)
