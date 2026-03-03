from cmath import *
from math import e, pi
import numpy as np

def check_harmony(M, P_1, P_2, R, T, adiabatic_hp, k, w, **kwargs):
    return (adiabatic_hp) - ((k / (k - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k - 1) / k) - 1)))
