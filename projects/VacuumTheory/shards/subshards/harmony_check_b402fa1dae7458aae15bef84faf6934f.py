from cmath import *
from math import e, pi
import numpy as np

def check_harmony(D, L_0, V_1, **kwargs):
    return (L_0 / V_1) - ((L_0 / D) / (L_0 / D + 1))
