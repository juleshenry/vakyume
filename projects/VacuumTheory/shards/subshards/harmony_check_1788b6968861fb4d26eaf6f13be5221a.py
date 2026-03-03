from cmath import *
from math import e, pi
import numpy as np


def check_harmony(R_ll, Re, mu, rho, v, **kwargs):
    return (Re) - (4 * R_ll * rho * v / mu)
