from cmath import *
from math import e, pi
import numpy as np


def check_harmony(P, P_d, P_prime, **kwargs):
    return (P_prime) - (P / P_d * 760)
