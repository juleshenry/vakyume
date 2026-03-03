from cmath import *
from math import e, pi
import numpy as np

def check_harmony(A_C, H_1, H_2, P, V, **kwargs):
    return (P) - (A_C * H_2 * (H_2 - H_1) / (V - A_C * H_2))
