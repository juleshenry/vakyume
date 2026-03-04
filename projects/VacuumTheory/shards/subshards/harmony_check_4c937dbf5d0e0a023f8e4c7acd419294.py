from cmath import *
from math import e, pi
import numpy as np


def check_harmony(H_2, KAPPA_1, P, **kwargs):
    return (P) - (KAPPA_1 * H_2**2)
