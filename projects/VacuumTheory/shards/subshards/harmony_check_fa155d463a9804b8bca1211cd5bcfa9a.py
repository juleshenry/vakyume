from cmath import *
from math import e, pi
import numpy as np


def check_harmony(K_1, K_2, alpha_1_2, **kwargs):
    return (alpha_1_2) - (K_1 / K_2)
