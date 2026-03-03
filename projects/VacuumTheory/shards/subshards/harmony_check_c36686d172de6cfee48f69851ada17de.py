from cmath import *
from math import e, pi
import numpy as np


def check_harmony(T, k, m, v, **kwargs):
    return (0.5 * m * v**2) - (1.5 * k * T)
