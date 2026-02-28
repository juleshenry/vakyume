from math import *
import numpy as np

def check_harmony(Re, h, mu, rho, v, w, **kwargs):
    return (Re) - ((2 * w * h * rho * v) / ((w + h) * mu))
