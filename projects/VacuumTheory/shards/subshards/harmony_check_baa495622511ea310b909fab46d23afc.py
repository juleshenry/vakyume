from cmath import *
from math import e, pi
import numpy as np

def check_harmony(P_1, P_2, adiabatic_power_watts, f, **kwargs):
    return (adiabatic_power_watts) - (f / 12 * ((P_2 / P_1) ** 0.286 - 1))
