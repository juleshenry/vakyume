from cmath import *
from math import e, pi
import numpy as np

def check_harmony(Suc_Pres, delta_P, oper_press, **kwargs):
    return (Suc_Pres) - (oper_press - delta_P)
