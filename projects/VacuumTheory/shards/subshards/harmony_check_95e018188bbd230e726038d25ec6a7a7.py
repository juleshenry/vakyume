from cmath import *
from math import e, pi
import numpy as np

def check_harmony(P, S_Th, S_p, T_e, T_i, p_c, p_s, **kwargs):
    return (S_p) - (S_Th * (P - p_s)*(460 + T_i) / ((P - p_c)*(460 + T_e) ))
