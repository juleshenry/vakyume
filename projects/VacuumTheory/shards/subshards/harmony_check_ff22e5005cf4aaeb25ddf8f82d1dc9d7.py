from cmath import *
from math import e, pi
import numpy as np

def check_harmony(C, S_1, S_2, **kwargs):
    return (S_1 ** -1) - (S_2 ** -1 + 1 / C)
