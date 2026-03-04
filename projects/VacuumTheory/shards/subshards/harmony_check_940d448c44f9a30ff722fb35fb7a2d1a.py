from cmath import *
from math import e, pi
import numpy as np


def check_harmony(P, S_0, S_p, T_e, T_i, p_0, p_c, p_s, **kwargs):
    return (S_0) - (
        S_p
        * ((P - p_0) * (460 + T_i) * (P - p_c) / (P * (P - p_s) * (460 + T_e))) ** 0.6
    )
