from cmath import *
from math import e, pi
import numpy as np

def check_harmony(P_s, Q_throughput, S_p, **kwargs):
    return (Q_throughput) - (S_p * P_s)
