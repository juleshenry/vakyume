from math import *
import numpy as np

def check_harmony(D, Re, mu, rho, v, **kwargs):
    return (Re) - (rho * D * v / mu)
