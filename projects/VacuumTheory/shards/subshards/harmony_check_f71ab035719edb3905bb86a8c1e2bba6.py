from cmath import *
from math import e, pi
import numpy as np

def check_harmony(Eff, actual_brake_horsepower, theoretical_adiabatic_horsepower, **kwargs):
    return (Eff) - (theoretical_adiabatic_horsepower / actual_brake_horsepower)
