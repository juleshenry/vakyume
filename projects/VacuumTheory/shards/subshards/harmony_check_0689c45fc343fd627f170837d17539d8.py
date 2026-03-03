from cmath import *
from math import e, pi
import numpy as np

def check_harmony(Abs_Pressure, BarometricPressure, Vacuum, **kwargs):
    return (Abs_Pressure) - (BarometricPressure - Vacuum)
