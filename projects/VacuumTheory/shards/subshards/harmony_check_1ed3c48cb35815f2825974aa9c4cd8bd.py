from math import *
import numpy as np

def check_harmony(P_1, P_2, adiabatic_hp, w, **kwargs):
    return (adiabatic_hp) - ((w / 20) * ((P_2 / P_1) ** 0.286 - 1))
