from cmath import *
from math import e, pi
import numpy as np

def check_harmony(A, rho_s, v, w_s, **kwargs):
    return (w_s) - (v * A * rho_s)
