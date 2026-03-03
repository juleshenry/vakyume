from cmath import *
from math import e, pi
import numpy as np

def check_harmony(V, r_h, t_h, w_h, **kwargs):
    return (w_h) - (r_h * V / t_h)
