from cmath import *
from math import e, pi
import numpy as np


def check_harmony(L_0, R, V_1, **kwargs):
    return (L_0 / V_1) - (R / (R + 1))
