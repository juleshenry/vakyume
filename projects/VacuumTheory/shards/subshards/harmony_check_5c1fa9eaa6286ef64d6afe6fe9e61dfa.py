from cmath import *
from math import e, pi
import numpy as np

def check_harmony(P, S_0, S_Th, p_0, **kwargs):
    return (S_Th) - (S_0 * (P / (P - p_0)) ** 0.6)
