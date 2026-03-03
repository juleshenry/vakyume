from cmath import *
from math import e, pi
import numpy as np


def check_harmony(C, D, L, P_p, mu, **kwargs):
    return (C) - (3.141592653589793 * D**4 / (128 * mu * L) * P_p)
