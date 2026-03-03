from cmath import *
from math import e, pi
import numpy as np


def check_harmony(P_s, V, t_e, w_j, **kwargs):
    return (t_e) - ((2.3 - 0.003 * P_s) * V / w_j)
