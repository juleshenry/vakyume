from cmath import *
from math import e, pi
import numpy as np


def check_harmony(Q_gas, SP_1, SP_2, S_p, V, t, **kwargs):
    return (t) - (V / S_p * log((SP_1 - Q_gas) / (SP_2 - Q_gas)))
