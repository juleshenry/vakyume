from cmath import *
from math import e, pi
import numpy as np


def check_harmony(Q_v, T_1, T_2, T_R, c_p, w_1, w_2, **kwargs):
    return (w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R)) - (12000 * Q_v)
