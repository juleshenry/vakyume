from math import *
import numpy as np

def check_harmony(Re, f, **kwargs):
    return (f) - (0.316 / Re ** (0.25))
