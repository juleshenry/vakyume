from math import *
import numpy as np

def check_harmony(Abs_Pressure, BarometricPressure, Vacuum, **kwargs):
    return (Abs_Pressure) - (BarometricPressure - Vacuum)
