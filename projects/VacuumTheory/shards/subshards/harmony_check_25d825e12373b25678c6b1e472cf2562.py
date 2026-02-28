from math import *
import numpy as np

def check_harmony(bhp, bhp_0, mu, rho, **kwargs):
    return (bhp) - (bhp_0 * (0.5 + 0.0155 * rho ** 0.84 * mu ** 0.16))
