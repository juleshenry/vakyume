from cmath import *
from math import e, pi
import numpy as np

def check_harmony(Eff, N_ES, N_t, T, **kwargs):
    return (N_t) - (N_ES / Eff ** T)
