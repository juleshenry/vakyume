from cmath import *
from math import e, pi
import numpy as np


def check_harmony(
    C_1, C_2, T_1, T_2, c_p, delta_h_c, delta_h_v, delta_t, m_b, w_v, **kwargs
):
    return (w_v) - (
        (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))
        / (delta_t * delta_h_v)
    )
