from cmath import *
from math import e, pi
import numpy as np

def check_harmony(B, L_N, V_0, **kwargs):
    return (L_N / V_0) - ((V_0 + B) / V_0)
