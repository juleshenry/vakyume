from math import *
import numpy as np

def check_harmony(bhp, c_p, delta_T, delta_h_i, f_a, rho, w_i, **kwargs):
    return (delta_T) - ((2545 * bhp + w_i * delta_h_i) / ( 8.02 * f_a * rho * c_p ))
