from cmath import *
from math import e, pi
import numpy as np


def check_harmony(H_1, H_2, KAPPA_2, P, **kwargs):
    return (P) - (KAPPA_2 * (H_2 - H_1))
