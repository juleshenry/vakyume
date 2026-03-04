from cmath import *
from math import e, pi
import numpy as np


def check_harmony(D, L, delta_P, mu, q, **kwargs):
    return (q) - (3.141592653589793 * (D**4) * delta_P / (128 * L * mu))
