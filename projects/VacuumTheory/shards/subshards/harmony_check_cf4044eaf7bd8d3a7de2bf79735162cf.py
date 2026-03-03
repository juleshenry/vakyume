from cmath import *
from math import e, pi
import numpy as np

def check_harmony(P_0_v, P_D, p_g, p_v_max, **kwargs):
    return (p_v_max) - (P_0_v * p_g / (P_D - P_0_v))
