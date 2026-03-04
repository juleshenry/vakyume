from cmath import *
from math import e, pi
import numpy as np


def check_harmony(M, P_c, T_c, mu_c, **kwargs):
    return (mu_c) - ((7.7 * (M**0.5) * P_c ** (2 / 3)) / T_c ** (1 / 6))
