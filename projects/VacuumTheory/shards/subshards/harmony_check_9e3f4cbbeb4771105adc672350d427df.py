from cmath import *
from math import e, pi
import numpy as np

def check_harmony(P, p_i, y_i, **kwargs):
    return (y_i) - (p_i / P)
