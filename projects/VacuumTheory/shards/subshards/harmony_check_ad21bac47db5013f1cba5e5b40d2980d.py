from cmath import *
from math import e, pi
import numpy as np

def check_harmony(L_c_P, Q_condensor_heat_duty, del_T, **kwargs):
    return (L_c_P) - (Q_condensor_heat_duty / (500 * del_T))
