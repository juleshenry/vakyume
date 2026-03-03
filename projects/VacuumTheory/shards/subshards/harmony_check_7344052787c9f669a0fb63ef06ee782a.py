from cmath import *
from math import e, pi
import numpy as np


def check_harmony(T_1, T_2, T_R, c_p, del_h_v, w_1, w_2, w_v, **kwargs):
    return (w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R)) - (w_v * del_h_v)
