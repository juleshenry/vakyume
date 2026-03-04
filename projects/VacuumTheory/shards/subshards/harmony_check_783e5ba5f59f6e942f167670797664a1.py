from cmath import *
from math import e, pi
import numpy as np


def check_harmony(F_s, t, t_c, **kwargs):
    return (t) - (t_c * F_s)
