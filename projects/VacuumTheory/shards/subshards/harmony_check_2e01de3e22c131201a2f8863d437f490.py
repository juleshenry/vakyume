from cmath import *
from math import e, pi
import numpy as np

def check_harmony(R_ll, h, w, **kwargs):
    return (R_ll) - (w * h / (2 * (w + h)))
