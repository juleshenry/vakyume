from cmath import *
from math import e, pi
import numpy as np

def check_harmony(P_1, P_2, T_1, T_2, V_1, V_2, **kwargs):
    return (P_1 * V_1 / T_1) - (P_2 * V_2 / T_2)
