from cmath import *
from math import e, pi
import numpy as np

def check_harmony(E_j, E_m, e, r, s, **kwargs):
    return (r) - (2.93 * (E_j * e) / (E_m * s))
