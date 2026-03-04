from cmath import *
from math import e, pi
import numpy as np


def check_harmony(C, C_1, C_2, D, L, P_p, mu, **kwargs):
    return (C) - (C_1 * (D**4 / (mu * L)) * P_p + C_2 * (D**3 / L))
