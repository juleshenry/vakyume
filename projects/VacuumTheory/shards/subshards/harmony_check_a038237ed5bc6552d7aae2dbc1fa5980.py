from cmath import *
from math import e, pi
import numpy as np


def check_harmony(P_c, p, p_i, p_nc, **kwargs):
    return (p_i / p_nc) - (p_i / (p - P_c))
