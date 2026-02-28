from math import *
import numpy as np

def check_harmony(_beta, mu, vel_grad, **kwargs):
    return (_beta) - (mu * vel_grad)
