from cmath import *
from math import e, pi
import numpy as np


def check_harmony(G, G_C, H, P, rho, **kwargs):
    return (P) - (G / (G_C * rho * H))
