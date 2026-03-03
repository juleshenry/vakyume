from cmath import *
from math import e, pi
import numpy as np


def check_harmony(P_0_1, P_0_2, alpha_12, **kwargs):
    return (alpha_12) - (P_0_1 / P_0_2)
