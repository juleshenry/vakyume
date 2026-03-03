from cmath import *
from math import e, pi
import numpy as np


def check_harmony(L, d, delta_P, f, q, rho, **kwargs):
    return (delta_P) - (2.15 * rho * f * L * q**2 / (d**5))
