from cmath import *
from math import e, pi
import numpy as np

def check_harmony(C, S_p, S_pump_speed, **kwargs):
    return (S_pump_speed) - ((S_p * C) / (S_p + C))
