from math import *
import numpy as np

def check_harmony(lambd, mu, rho, v_a, **kwargs):
    return (mu) - (0.35 * rho * lambd * v_a)
