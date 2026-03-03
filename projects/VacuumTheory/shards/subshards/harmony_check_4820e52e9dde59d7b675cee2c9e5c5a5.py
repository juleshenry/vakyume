from cmath import *
from math import e, pi
import numpy as np


def check_harmony(L, d, delta_P, mu, v, **kwargs):
    return (delta_P) - (0.0345 * mu * L * v / d**2)
