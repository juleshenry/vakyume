from cmath import *
from math import e, pi
import numpy as np


def check_harmony(P_0_V, P_D, P_v_0, S_B, S_D, p_b, p_g, p_v_max, **kwargs):
    return (p_v_max) - (
        S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
    )
