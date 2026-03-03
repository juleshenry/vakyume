from cmath import *
from math import e, pi
import numpy as np

def check_harmony(F_ext, M, a_CM, **kwargs):
    return (F_ext) - (M * a_CM)
