from math import *
import numpy as np

def check_harmony(M, P, R, T, rho, **kwargs):
    return (rho) - (P * M / (R * T))
