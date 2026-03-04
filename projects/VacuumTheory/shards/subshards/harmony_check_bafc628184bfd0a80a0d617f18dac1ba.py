from cmath import *
from math import e, pi
import numpy as np


def check_harmony(A, dV_dt, delta_P, m, mu, r, r_M, **kwargs):
    return (dV_dt) - ((A * delta_P) / (mu * (m / A) * r * delta_P + r_M))
