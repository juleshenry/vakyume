from cmath import *
from math import e, pi
import numpy as np

def check_harmony(P_1, P_2, S_p, V, t, **kwargs):
    return (t) - (V / S_p * log(P_1 / P_2))
