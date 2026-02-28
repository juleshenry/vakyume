from math import *
import numpy as np

def check_harmony(L, d, delta_P, f, g, rho, v, **kwargs):
    return (delta_P) - (4.31 * rho * f * L * v ** 2 / (2 * d * g))
