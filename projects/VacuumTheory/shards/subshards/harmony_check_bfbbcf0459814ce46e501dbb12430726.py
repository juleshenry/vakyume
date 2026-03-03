from cmath import *
from math import e, pi
import numpy as np

def check_harmony(M_1, M_2, P_0_1, P_0_2, a_M_12, **kwargs):
    return (a_M_12) - ((P_0_1) / (P_0_2) * (M_2 / M_1) ** 0.4)
