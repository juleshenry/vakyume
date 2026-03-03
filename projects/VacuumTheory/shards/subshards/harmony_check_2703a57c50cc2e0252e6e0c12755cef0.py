from cmath import *
from math import e, pi
import numpy as np

def check_harmony(K_1, K_2, x_1, x_2, y_1, y_2, **kwargs):
    return (K_1 / K_2) - (y_1 * x_2 / (y_2 * x_1))
