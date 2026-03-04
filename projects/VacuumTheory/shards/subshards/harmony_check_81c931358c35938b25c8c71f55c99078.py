from cmath import *
from math import e, pi
import numpy as np


def check_harmony(P_c, n_i, n_nc, p, p_i, **kwargs):
    return (n_i / n_nc) - (p_i / (p - P_c))
