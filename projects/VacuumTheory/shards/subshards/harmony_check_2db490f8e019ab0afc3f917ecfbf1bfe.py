from cmath import *
from math import e, pi
import numpy as np

def check_harmony(C, P_1, P_2, Q_throughput, **kwargs):
    return (C) - (Q_throughput / (P_1 - P_2))
