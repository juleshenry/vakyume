from cmath import *
from math import e, pi
import numpy as np

def check_harmony(delta, lambd, psi, **kwargs):
    return (lambd) - (3.141592653589793 * delta ** 2 * psi * 2 ** 0.5)
