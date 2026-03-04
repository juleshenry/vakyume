from cmath import *
from math import e, pi
import numpy as np


def check_harmony(D_r, sig_R, w, **kwargs):
    return (sig_R) - (0.00436 * D_r * w)
