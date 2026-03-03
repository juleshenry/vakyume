from cmath import *
from math import e, pi
import numpy as np


def check_harmony(A, dV_dt, delta_P, mu, r_c, s, tau, **kwargs):
    return (dV_dt) - ((A * delta_P ** (1 - s)) / (mu * tau * r_c))
